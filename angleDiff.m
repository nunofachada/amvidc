%
% angleDiff function - determine smallest angle between two directions
%
% Parameters:
%       v1 - n x 1, vector representing first direction
%       v2 - n x 1, vector representing second direction
% Output:
%    delta - angle in radians between first and second directions
%
function delta = angleDiff(v1, v2)

% Determine number of dimensions
numDims = size(v1, 1);

% Make sure vectors are in 1st or 2nd quadrants (from a 2D perspective,
% although this should work for m-dimensions)
if  v1(numDims, 1) < 0
    v1 = -1 * v1;
end;
if  v2(numDims, 1) < 0
    v2 = -1 * v2;
end;

% Obtain angle between vectors (thus, between the directions they
% represent)
cosDelta = dot(v1, v2) / (norm(v1) * norm(v2));
delta = acos(cosDelta);
%delta = 1 - cosDelta;

end
