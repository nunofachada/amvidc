function idx = initClust(data, minClustSize, distance, numex)
% INITCLUST Performs very simple initial clustering based on AHC with 
% single linkage (nearest neighbor) and user defined distance. Each sample 
% is associated with the same cluster of its nearest point.
%
% idx = INITCLUST(data, minClustSize, distance, numex)
%
% Parameters:
%            data - data to cluster
%    minClustSize - minimum size of each cluster
%        distance - type of distance to consider (as supported my matlab
%                   pdist function), default is 'seuclidean'
%           numex - number of clusters which are allowed to have less than
%                   the minimum size (default = 1)
% Output:
%             idx - clustering result
%

%  N. Fachada
%  Instituto Superior Técnico, Lisboa, Portugal

% Default distance
if nargin < 3, distance = 'seuclidean'; end;
% Default number of exceptions
if nargin < 4, numex = 1; end;
    
% Determine number of samples
numSamples = size(data, 1);

% Obtain distance matrix
distMatrix = squareform(pdist(data, distance));

% Make zeros in distance matrix infinite
distMatrix = distMatrix + diag(ones(1, size(distMatrix, 1)) * inf);

% Initial clustering is one sample per cluster
idx = (1:numSamples)';

% Connect clusters until all (minus numex clusters) are connected
while sum(histc(idx, unique(idx)) < minClustSize) > numex
    for i=1:numSamples
        % Check if sample is already connected x times
        timesConnected = max(size(find(distMatrix(i, :) == inf))) - 1;
        %fprintf('%d connected %d times: ', i, timesConnected);
        if timesConnected >= minClustSize - 1
            %fprintf('%d already connected enough times!\n', i);
            continue;
        end;
        % If not, connect to closer sample
        [~, index] = min(distMatrix(i, :));
        distMatrix(i, index) = inf;
        distMatrix(index, i) = inf;
        %fprintf('%d connected to %d\n', i, index);
        oldClust = idx(index);
        idx(idx == oldClust) = idx(i);
    end;
end;

% Return idx, normalized
idx = idxNormalize(idx);