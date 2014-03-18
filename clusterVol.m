% 
% clusterVol function - Determine the volume of the cluster formed by the
% given set of points. The volume is determined using the convex hull
% formed by the cluster, of the minimum volume ellipsoid (mve) formed by 
% the cluster.
% 
% Parameters:
%       points - m x n, with m samples and n dimensions
%         type - 'ellipsoid' or 'convhull'
% zeroVolValue - value to assign to volume if given points are not enough
%                to calculate a volume
%          tol - tolerance for ellipsoid volume
% Output:
%  volCluster - Volume of the cluster formed by the given set of points.
%
function volCluster = clusterVol(points, type, zeroVolValue, tol)

% How many points are in cluster
sizeCluster = size(points, 1);
% How many dimensions are at stake
numDims = size(points, 2);
% Add path to external functions
%addpath('external/');

% Check if there are enough points to calculate a volume
if sizeCluster < numDims + 1
    volCluster = zeroVolValue;
    return;
end;

% Determine what type of volume to use
if strcmp(type, 'ellipsoid')
    % Ellipsoid, use MinVolEllipse from Nima Moshtagh
    [A, ~] = mve.MinVolEllipse(points', tol);
    volCluster = det(inv(A));
elseif strcmp(type, 'convhull')
    % Convex hull, use matlab native function
    [~, volCluster] = convhulln(points,{'QJ','Pp'});
else
    error('Unknown type of volume.');
end;

