function wd = weightedDirection(vectors, weights)
% WEIGHTEDDIRECTION Return a weighted direction given several vectors and 
% correspondent weights.
%
% wd = WEIGHTEDDIRECTION(vectors, weights)
%
% Parameters:
%     vectors - m x n, m dimensions, n vectors (normalized)
%     weights - 1 x n, weight of each vector
% Output:
%          wd - m x 1, normalized weighted direction
%

%  N. Fachada
%  Instituto Superior Técnico, Lisboa, Portugal

% Determine number of dimensions
numDims = size(vectors, 1);

% Determine number of vectors
numVectors = size(vectors, 2);

% Make sure vectors are in 1st or 2nd quadrants (from a 2D perspective,
% although this should work for m-dimensions)
for i=1:numVectors
    if  vectors(numDims, i) < 0
        vectors(:, i) = -1 * vectors(:, i);
    end;
end;

% Determine weighted global direction
wd = vectors * weights';
wd = wd / norm(wd);
