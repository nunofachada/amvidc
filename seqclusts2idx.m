function idx = seqclusts2idx(seqclusts)
% SEQCLUSTS2IDX Convert a sequence of number of clusters to a 
% sample to cluster map (an idx, as returned by kmeans, etc).
%
% idx = SEQCLUSTS2IDX(seqclusts)
%
% Parameters:
%   seqclusters - p x 1, number of samples per cluster
% Output:
%           idx - m x 1, maps each sample with a cluster
%

%  N. Fachada
%  Instituto Superior Técnico, Lisboa, Portugal

% Number of samples
numSamples = sum(seqclusts);
% Number of clusters
numClusts = max(size(seqclusts));
% Allocate space to idx
idx = zeros(numSamples, 1);
% Helper variable
idxIndexes = cumsum(seqclusts);

for i=1:numClusts
    if i==1
        from = 1;
    else
        from = idxIndexes(i - 1) + 1;
    end;
    to = idxIndexes(i);
    idx(from:to, 1) = i;
end;
