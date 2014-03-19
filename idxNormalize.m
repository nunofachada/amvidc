function normalizedIdx = idxNormalize(idx)
% IDXNORMALIZE Reshifts idx cluster id's to be within 1 and the total 
% number of clusters.
%
% normalizedIdx = IDXNORMALIZE(idx)
%
% Parameters:
%            idx - typical idx clustering result
% Output:
%  normalizedIdx - normalized idx clustering result
%

%  N. Fachada
%  Instituto Superior Técnico, Lisboa, Portugal

idxUniques = [unique(idx) zeros(size(unique(idx), 1) ,1)];
normalizedIdx = zeros(size(idx, 1), 1);

clustCount = 0;
for i=1:size(idx, 1)
    index = find(idx(i) == idxUniques(:, 1));
    if idxUniques(index, 2) == 0
        clustCount = clustCount + 1;
        idxUniques(index, 2) = clustCount;
    end;
    normalizedIdx(i, 1) = idxUniques(index, 2);
end;

