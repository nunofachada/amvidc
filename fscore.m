%
% fscore function - Evaluate clustering result by comparing it with 
% correct result. Assumes that correct result has observations of same 
% class in sequence.
%
% Parameters:
%              idx - Clustering result (e.g. such as given by kmeans)
%       numclasses - Number of clusters/classes
%  numclassmembers - Vector with size of each cluster (or a scalar if all
%                    clusters are of the same size)
% Output:
%             eval - clustering fscore of idx for given classes
%
function eval = fscore(idx, numclasses, numclassmembers)

% If numclassmembers is a scalar, transform it in a vector
if max(size(numclassmembers)) == 1
    numclassmembers = numclassmembers * ones(1, numclasses);
end;

% Set each class fscore to zero
fscores = zeros(1, numclasses);

% Determine the size of each cluster
numclustmembers = histc(idx, 1:numclasses);
    
% Determine F-Score for each class
for i=1:numclasses
    
    % Determine starting index of idx for current class
    if i > 1
        indexFrom = sum(numclassmembers(1:(i-1))) + 1;
    else
        indexFrom = 1;
    end;
    % Determine endind index of idx for current class
    indexTo = sum(numclassmembers(1:i));

    % Must experiment for each cluster
    classPointsInClusters = histc(idx(indexFrom:indexTo), 1:numclasses);
    for j=1:numclasses
        recall = classPointsInClusters(j) / numclassmembers(i);
        precision = classPointsInClusters(j) / numclustmembers(j);
        fscore_class = 2 * recall * precision / (recall + precision);
        % If this is the highest fscore for this class, retain it
        if (fscore_class > fscores(i))
            fscores(i) = fscore_class;
        end;
    end;
end;

% Determine final fscore
eval = sum(fscores .* numclassmembers) / max(size(idx));

