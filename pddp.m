% 
% pddp function - Perform PDDP (principal direction divisive clustering)
%                 on input data.
% 
% Parameters:
%       data -  m x n, with m samples and n dimensions
%  maxClusts - maximum number of clusters to form (optional, default is inf)
%  minPoints - minimum number of points in each cluster (optional, default
%              is ceil((size(data, 2)+1)/2) )
%      debug - Show debug info? (optional, default is false)
% Output:
%         idx - final clustering result
%
function idx = pddp(data, maxClusts, minPoints, debug)

% Check if maximum number of clusters is given as argument
if nargin < 2, maxClusts = inf; end;
% Check if minimum number of points in each cluster is given
if nargin < 3, minPoints = ceil((size(data, 2) + 1)/2); end;
% Check if debug info is given
if nargin < 4, debug = 0; end;

% Get number of samples and number of dimensions
numSamples = size(data, 1);

% Pre-allocate idx
idx = ones(numSamples, 1);

% Start algorithm
while true
    % Get largest cluster
    largest = mode(idx);
    % Get indexes of points in largest cluster
    iIndexes = find(idx == largest);
    % How many points does the largest cluster contain?
    numObs = size(iIndexes, 1);
    % Print some info
    debugf(sprintf('Largest cluster has %d observations, data has %d dimensions, minPoints is %d.\n', numObs, size(data, 2), minPoints), debug);
    % See if its big enough to make a volume
    if  numObs > minPoints * 2
        % Perform
        [~, SCORE] = princomp(data(iIndexes, :), 'econ');
        newCluster = max(idx) + 1;
        newIndexes = iIndexes(SCORE(:, 1) > 0);
        numObsNewClust = size(newIndexes, 1);
        debugf(sprintf(' - Number of observations to go to new cluster: %d\n', numObsNewClust), debug);
        if min(numObsNewClust, numObs - numObsNewClust) < minPoints
            debugf(' - One of the clusters is not big enough, terminating algorithm...\n', debug);
            break;
        end;
        idx(newIndexes) = newCluster;
        if newCluster == maxClusts
            debugf(' - Reached maximum number of clusters, terminating algorithm...\n', debug);
            break;
        end;
    else
        % Can't divide no more, end algorithm
        debugf(' - Largest cluster is not divisible, terminating algorithm...\n', debug);
        break;
    end;
end;

end % pddp

% Debug helper function
function debugf(str, debug)
    if debug
        fprintf(str);
    end;
end


    