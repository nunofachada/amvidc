function [idx, idx_all, vol_all, dir_all, diss_all, min_all] = ...
    clusterdata_amvidc(data, k, idx, varargin)
% CLUSTERDATA_AMVIDC Perform agglomerative hierarchical clustering (AHC) 
% using Minimum Volume Increase (MVI) and Minimum Direction Change (MDC) 
% criteria.
%
% [idx idx_all vol_all dir_all diss_all min_all] = 
%   CLUSTERDATA_AMVIDC(data, k, idx, varargin)
%
% Parameters:
%        data - m x n, with m samples and n dimensions
%           k - number of clusters
%         idx - m x 1, initial clustering
%    varargin - optional parameters as key-value pairs:
%             'volume' - 'ellipsoid', 'convhull' (default)
%                'tol' - for ellipsoid volume only (tolerance for minimum
%                        volume ellipse calculation), default = 0.01
%          'dirweight' - final direction weight (default is 0, i.e.,
%                        cluster directionallity is ignored)
%           'dirpower' - convergence power of effective dirweight > 0
%                        (higher powers make convergence steeper and
%                        ocurring more to the end).
%            'dirtype' - 'pca', 'svd' (default)
%                'nvi' - allow negative volume increase? default is true
%           'loglevel' - log level: 0 (show all) to 4 (only show critical
%                        errors), default is 3 (show warnings)
%
% Output:
%         idx - final clustering result
%     idx_all - clusters formed during each step of AHC
%     vol_all - volume changes calculated in each step of AHC
%     dir_all - direction changes calculated in each step of AHC
%    diss_all - dissimilarity matrix (vol + dir changes) in each step of AHC
%     min_all - minimum change (MVI+MDC) in each step of AHC
%

%  N. Fachada
%  Instituto Superior Técnico, Lisboa, Portugal

% Set global logging level structure
global loglevels;
global loglevel;
loglevels = struct('ALL', 0, 'DEBUG', 1, 'INFO', 2, 'WARN', 3, 'ERROR', 4);

% Parse optional arguments or get default values for respective parameters
[volume, tol, dirweight, dirpower, dirtype, nvi, loglevel] = ...
    parseOptionalArgs(varargin);

% Current clustering step
ind = 1;

% Obtain g (number of clusters)
g = size(unique(idx), 1);

% Initialize clustering history
numClusteringSteps = g - k;
idx_all = cell(1, numClusteringSteps);
dir_all = cell(1, numClusteringSteps);
diss_all = cell(1, numClusteringSteps);
min_all = zeros(1, numClusteringSteps);
vol_all = cell(1, numClusteringSteps);
idx_all{ind} = idx;


% % % %
% % % % Determine initial volumes, and save time later
% % % %

logmsg('Determining initial volumes....\n', loglevels.INFO);

% Init. volume variables
vol = inf(g, g);
max_vol = -inf; % This is used for normalization later

% Compute volume increase of all possible cluster unions
for i=1:g-1
    for j=i+1:g
        % Get indexes of samples in cluster i
        iIndexes = idx == i;
        % Get indexes of samples in cluster j
        jIndexes = idx == j;
        % Get samples in cluster i
        pointsInCluster1 = data(iIndexes, :);
        % Get samples in cluster j
        pointsInCluster2 = data(jIndexes, :);
        % Get samples in new cluster (i and j combined)
        pointsInNewCluster = [pointsInCluster1 ; pointsInCluster2];
        % Determine the volume difference between the new cluster and
        % the clusters i and j
        vol(i, j) = ...
            clusterVol(pointsInNewCluster, volume, inf, tol) - ...
            clusterVol(pointsInCluster1, volume, 0, tol) - ...
            clusterVol(pointsInCluster2, volume, 0, tol);
        % Check for negative volume increase
        if vol(i, j) < 0
            if ~nvi
                vol(i, j) = inf;
            end;
            logmsg(['Volume increase =< 0 detected in iteration ' int2str(ind)], loglevels.WARN);
        end;
        % Make sure volume dissimilarity matrix is simetric
        vol(j, i) = vol(i, j);
        % If this is the maximum volume, keep it (required for later
        % normalization)
        if vol(i, j) > max_vol && ~isinf(vol(i, j))
            max_vol = vol(i, j);
        end;
    end;
end;
% Normalize volume
vol_norm = vol / max_vol;
% Update volume history
vol_all{1} = vol;

% % % %
% % % % Clustering Stage using AHC with MVI+DC linkage
% % % %

logmsg('Entering cycle....\n', loglevels.DEBUG);

% Do this while current number of clusters g is higher than final number of
% clusters k
while g > k
    
    logmsg(sprintf('\n-----------\n *iterations left: #%d* \n-----------\n', g-k), loglevels.INFO);
    
    % Init. direction variables
    direction = inf(g, g);
    dir_norm = direction;
    max_direction = 0; % for normalization purposes later
    
    % Perform direction calculations only if dirweight > 0
    if dirweight > 0
        % Compute current global direction
        globalDirection = getPrincipalDirection(data, idx, g, dirtype);
        % Compute direction difference of new cluster to global direction
        for i=1:g-1
            for j=i+1:g
                % Get indexes of samples in cluster i
                iIndexes = idx == i;
                % Get indexes of samples in cluster j
                jIndexes = idx == j;
                % Get samples in cluster i
                pointsInCluster1 = data(iIndexes, :);
                % Get samples in cluster j
                pointsInCluster2 = data(jIndexes, :);
                % Get samples in new cluster (i and j combined)
                pointsInNewCluster = [pointsInCluster1 ; pointsInCluster2];
                % Determine the principal direction of new cluster
                newClustPD = getPrincipalDirection(pointsInNewCluster, ones(size(pointsInNewCluster, 1), 1), 1, dirtype);
                % Determine the diretion difference between the new cluster
                % and the global direction
                direction(i, j) = angleDiff(newClustPD, globalDirection);
                direction(j, i) = direction(i, j);
                % Keep maximum direction for normalization purposes
                if direction(i, j) > max_direction
                    max_direction = direction(i, j);
                end;
            end;
        end;
        % Normalize direction
        dir_norm = direction / max_direction;
    end;
    % Update direction history
    dir_all{ind} = direction;
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % Merge clusters which have minimum volume + direction change %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
    % Determine effective direction weight for current iteration
    dirWeightEffective = dirweight * (k / (g - 1))^dirpower;
    volWeightEffective = 1 - dirWeightEffective;
    
    % Determine dissimilarity matrix (MVI+DC) for current iteration
    dissMatrix = zeros(size(vol_norm));
    if dirWeightEffective > 0
        dissMatrix = dirWeightEffective * dir_norm;
    end;
    if volWeightEffective > 0
        dissMatrix = dissMatrix + volWeightEffective * vol_norm;
    end;
    
    % Keep dissimilarity matrix history
    diss_all{ind} = dissMatrix;
    
    % Determine minimum values for each column of the dissimilarity matrix
    % and the respective row indexes
    [mins minRows] = min(dissMatrix);
    
    % Find the dissimilarity matrix absolute minimum and the respective
    % column index
    [minimum minCol] = min(mins);
    
    % In case there are several equal minimums, chose the column index of
    % the first one
    minCol = minCol(1, 1);
    
    % Determine the row index of the chosen absolute minimum
    minRow = minRows(minCol);
    
    % The index of the new cluster will be the minimum value between the
    % absolute minimum column and row indexes
    newCluster = min(minCol, minRow);
    
    % The index of the old cluster will be the maximum value between the
    % absolute minimum column and row indexes
    oldCluster = max(minCol, minRow);
    
    % The samples associated with the old cluster will now be associated
    % with the new cluster
    idx(idx == oldCluster, 1) = newCluster;
    
    % The samples associated with the cluster holding the largest id (=g)
    % will be associated with the index of the old cluster
    idx(idx == g, 1) = oldCluster;
    
    % Update absolute minimum history
    min_all(ind) = minimum;
    
    % Print info
    logmsg(sprintf('| Num. clusts=%d | weighted mvi+mdc=%g  | mvi=%g | mdc=%g | dirweight=%f |.\n', max(size(unique(idx))), minimum, min(min(vol)), min(min(direction)), dirWeightEffective), loglevels.INFO);
    
    % Stop if the end is reached
    if g - 1 == k
        break;
    end;
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % Update volume dissimilarity matrix for next iteration %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
    % 1 - Invalidate new cluster row/column
    vol(newCluster, :) = inf;
    vol(:, newCluster) = inf;
    
    % 2 - Replace old cluster with last row/column
    if oldCluster ~= g
        % Replace by rows
        vol = [ ...
            vol(1:(oldCluster - 1), :) ; ...
            vol(g, :); ... % Higher id cluster goes to where old cluster was
            vol((oldCluster+1):(g-1), :) ...
            ];
        % Replace by columns
        vol = [ ...
            vol(:, 1:(oldCluster-1)) ...
            vol(:, g) ... % Higher id cluster goes to where old cluster was
            vol(:, (oldCluster+1):(g-1)) ...
            ];
        % Put inf in intersection of old cluster / last row/column
        vol(oldCluster, oldCluster) = inf;
    else
        % This case is simpler, old cluster add index = g, thus we only
        % need to eliminate last row and column
        vol = vol(1:(g - 1), :);
        vol = vol(:, 1:(g - 1));
    end;
    
    % 3 - Determine new volume differences for new cluster
    i = newCluster;
    % Get indexes of samples in cluster i
    iIndexes = idx == i;
    % Get samples in cluster i
    pointsInCluster1 = data(iIndexes, :);
    % Determine volume differences between new cluster and existing
    % clusters
    for j=1:(g-1)
        % If row and column are the same, value is inf
        if i==j
            vol(i, j) = inf;
            continue;
        end;
        % Get indexes of samples in cluster j
        jIndexes = idx == j;
        % Get samples in cluster j
        pointsInCluster2 = data(jIndexes, :);
        % Get samples in new cluster (clusters i and j combined)
        pointsInNewCluster = [pointsInCluster1 ; pointsInCluster2];
        % Determine the volume difference between the new cluster and the
        % clusters i and j
        vol(i, j) = ...
            clusterVol(pointsInNewCluster, volume, inf, tol) - ...
            clusterVol(pointsInCluster1, volume, 0, tol) - ...
            clusterVol(pointsInCluster2, volume, 0, tol);
        % Check for negative volume increase
        if vol(i, j) < 0
            if ~nvi
                vol(i, j) = inf;
            end;
            logmsg(['Volume increase =< 0 detected in iteration ' int2str(ind)], loglevels.WARN);
        end;
        % Make sure volume dissimilarity matrix is simetric
        vol(j, i) = vol(i, j);
        % If this is the maximum volume, keep it (required for later
        % normalization)
        if vol(i, j) > max_vol && ~isinf(vol(i, j))
            max_vol = vol(i, j);
        end;
    end;
    
    % Normalize volume
    vol_norm = vol / max_vol;
    % Update volume history
    vol_all{ind + 1} = vol;
    
    % One less cluster
    g = g - 1;
    ind = ind + 1;
    idx_all{ind} = idx;
    
end;

% Finish algorithm

end % main function

% % % % % % % % % % % % % % % % % % % % % % % %
% Helper function to parse optional arguments %
% % % % % % % % % % % % % % % % % % % % % % % %
function [volume, tol, dirweight, dirpower, dirtype, nvi, loglevel] = ...
    parseOptionalArgs(args)

global loglevels;

% Check if any arguments are given
numArgs = size(args, 2);
if numArgs > 0
    % arguments must come in pairs (maximum of 7 key-value pairs)
    if mod(numArgs, 2) == 0 && numArgs <= 14
        % Parse arguments
        for i=1:2:numArgs-1
            if strcmp(args{i}, 'volume')
                % Volume type
                volume = args{i + 1};
            elseif strcmp(args{i}, 'tol')
                % Tolerance for minimum volume ellipse calculation
                tol = args{i + 1};
            elseif strcmp(args{i}, 'dirweight')
                % Final direction weight
                dirweight = args{i + 1};
            elseif strcmp(args{i}, 'dirpower')
                % Direction power
                dirpower = args{i + 1};
            elseif strcmp(args{i}, 'dirtype')
                % Direction type
                dirtype = args{i + 1};
            elseif strcmp(args{i}, 'nvi')
                % Allow zero or negative volume increase
                nvi = args{i + 1};
            elseif strcmp(args{i}, 'loglevel')
                % Debug mode
                loglevel = args{i + 1};
            else
                % Oops... unknown parameter
                error('Unknown parameter "%s"', args{i});
            end;
        end;
    else
        % arguments must come in pairs (maximum of 7 key-value pairs)
        error('Incorrect number of optional arguments.');
    end;
end;

% Check if parameters are set, if not, use defaults
if ~exist('volume', 'var')
    % Volume type
    volume = 'convhull';
end;
if ~exist('tol', 'var')
    % Tolerance for minimum volume ellipse calculation
    tol = 0.01;
end;
if ~exist('dirweight', 'var')
    % Final direction weight
    dirweight = 0;
end;
if ~exist('dirpower', 'var')
    % Direction power
    dirpower = 2;
end;
if ~exist('dirtype', 'var')
    % Direction type
    dirtype = 'svd';
end;
if ~exist('nvi', 'var')
    % Allow zero or negative volume increase
    nvi = true;
end;
if ~exist('loglevel', 'var')
    % Distance limit
    loglevel = loglevels.WARN;
end;

end % parseOptionalArgs

% % % % % % % % % % % % % % % % % % % % % % % %
% Logger helper function                      %
% % % % % % % % % % % % % % % % % % % % % % % %

function logmsg(msg, level)

global loglevels;
global loglevel;

if level >=  loglevel
    if level < loglevels.WARN
        fprintf(msg);
    else
        if level == loglevels.WARN
            warning(msg);
        else
            error(msg);
        end;
    end
end

end % logmsg
