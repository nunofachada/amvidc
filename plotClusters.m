function handle = plotClusters(data, dims, idx_marker, idx_encircle, groupType, handle)
% PLOTCLUSTERS Plot a set of clusters with different markers and with 
% different encircle areas for points of different clusters. This function 
% accepts idxs returned by kmeans function and such, or an array of size 
% equal to the number of clusters, where each value is the number of points 
% in each cluster (thus, data will be ordered by clusters).
%
% handle = PLOTCLUSTERS(data, dims, idx_marker, idx_encircle, groupType, handle)
%
% Parameters:
%           data - m x n, with m samples and n dimensions
%           dims - number of dimensions to plot (2 or 3)
%     idx_marker - cluster indexes for markers grouping
%   idx_encircle - cluster index for encircle grouping (default = idx_marker)
%      groupType - group clusters using 'convhull' (default if idx_encircle
%                  is given), 'ellipsoid' or 'none' (default if idx_encircle is not
%                  given)
%         handle - figure handle (default = new handle)
% Output:
%      handle - plot figure handle
%

%  N. Fachada
%  Instituto Superior Técnico, Lisboa, Portugal

% Parse optional parameters
if nargin < 4, idx_encircle = idx_marker; end;
if nargin < 5
    if nargin < 4, groupType = 'none'; else groupType = 'convhull'; end;
end;
if nargin < 6, handle = figure(); end;
    
% Cut data into specified dimensions, if necessary
data = data(:, 1:dims);

% Markers
markers = [{'+b'} {'om'} {'*r'} {'xr'} {'sg'} {'ok'} {'*k'} {'pg'} {'+m'} {'oc'} {'^b'} {'hc'} {'sb'} {'vm'} {'<r'} {'hr'} {'>g'} {'pk'} {'*k'} {'.g'} {'py'} {'>c'} {'<b'} {'.c'} {'.b'} {'hm'} {'>r'} {'.r'} {'+g'} {'dk'} {'.k'} {'dg'} {'sy'} {'pc'} {'vb'} {'^c'} {'+b'} {'om'} {'*r'} {'xr'} {'sg'} {'ok'} {'*k'} {'pg'} {'+m'} {'oc'} {'^b'} {'hc'} {'sb'} {'vm'} {'<r'} {'hr'} {'>g'} {'pk'} {'*k'} {'.g'} {'py'} {'>c'} {'<b'} {'.c'} {'.b'} {'hm'} {'>r'} {'.r'} {'+g'} {'dk'} {'.k'} {'dg'} {'sy'} {'pc'} {'vb'} {'^c'} {'+b'} {'om'} {'*r'} {'xr'} {'sg'} {'ok'} {'*k'} {'pg'} {'+m'} {'oc'} {'^b'} {'hc'} {'sb'} {'vm'} {'<r'} {'hr'} {'>g'} {'pk'} {'*k'} {'.g'} {'py'} {'>c'} {'<b'} {'.c'} {'.b'} {'hm'} {'>r'} {'.r'} {'+g'} {'dk'} {'.k'} {'dg'} {'sy'} {'pc'} {'vb'} {'^c'} {'+b'} {'om'} {'*r'} {'xr'} {'sg'} {'ok'} {'*k'} {'pg'} {'+m'} {'oc'} {'^b'} {'hc'} {'sb'} {'vm'} {'<r'} {'hr'} {'>g'} {'pk'} {'*k'} {'.g'} {'py'} {'>c'} {'<b'} {'.c'} {'.b'} {'hm'} {'>r'} {'.r'} {'+g'} {'dk'} {'.k'} {'dg'} {'sy'} {'pc'} {'vb'} {'^c'}];
markers = [markers markers markers markers markers markers];

% Prepare figure
clf;
hold on;
grid on;

% Get number of clusters to be drawn using markers
numClusters = getClusterData(data, idx_marker);

% Plot clusters with markers
for i=1:numClusters
    % Get samples in cluster i
    clusterSamples = getClusterData(data, idx_marker, i);
    % Plot cluster
    if dims == 3
        % In 3D
        plot3(clusterSamples(:, 1), clusterSamples(:, 2), clusterSamples(:, 3), markers{i});
    else
        % In 2D
        plot(clusterSamples(:, 1), clusterSamples(:, 2), markers{i});
    end;
end;

% If cluster encirclement is not disabled, perform encirclement
if ~strcmp(groupType, 'none')
    
    % Get number of clusters to be encircled
    numClusters = getClusterData(data, idx_encircle);
    
    % Plot clusters encirclements
    for i=1:numClusters
        
        % Get samples in cluster i
        clusterSamples = getClusterData(data, idx_encircle, i);
        
        if size(clusterSamples, 1) == 2
            % If there are only two samples, draw a line connecting them
            if dims == 2
                line( ...
                    [clusterSamples(1, 1) clusterSamples(2, 1)], ...
                    [clusterSamples(1, 2) clusterSamples(2, 2)], ...
                    'Color', 'k' ...
                );
            elseif dims == 3
                line( ...
                    [clusterSamples(1, 1) clusterSamples(2, 1)], ...
                    [clusterSamples(1, 2) clusterSamples(2, 2)], ...
                    [clusterSamples(1, 3) clusterSamples(2, 3)], ...
                    'Color', 'k' ...
                );
            end;
        elseif size(clusterSamples, 1) > 2
            % If there are more than two samples, draw encirclement
            if strcmp(groupType, 'convhull')
                % Convex hull encirclement
                clear k;
                k = convhull(clusterSamples, 'simplify', true);
                if dims == 2
                    plot(clusterSamples(k, 1), clusterSamples(k, 2), 'k-');
                elseif dims == 3
                    trimesh(k, ...
                        clusterSamples(:,1), ...
                        clusterSamples(:,2), ...
                        clusterSamples(:,3), ...
                        'EdgeColor', 'k', 'FaceAlpha', 0.5 ...
                     );
                end;
            elseif strcmp(groupType, 'ellipsoid')
                % Ellipsoid encirclement
                [A , c] = MinVolEllipse(clusterSamples', 0.01);
                Ellipse_plot(A, c);
                %grid on;
            end;
        end;
    end;
end;

% Add labels
xlabel('x');
ylabel('y');
if dims == 3
    zlabel('z');
end;

return;

% Helper function, which identifies how idx is organized, and
% returns number of clusters (if no cluster is specified) OR
% points in specified cluster.
function clusterData = getClusterData(data, idx, varargin)
% Determine what idx contains
if max(size(idx)) == size(data, 1)
    % % % idx contains cluster to which each sample belongs to
    % Determine number of clusters
    clusterTags = unique(idx);
    numClusters = max(size(unique(clusterTags)));
    % If varargin is not specified, return number of clusters
    if size(varargin, 2) == 0
        clusterData = numClusters;
        return;
    end;
    % Otherwise return points in cluster specified in varargin
    clusterId = varargin{1};
    clusterData = data(idx == clusterTags(clusterId), :);
else
    % % % idx is array with size equal to the number of clusters
    % Determine number of clusters
    numClusters = max(size(idx));
    % If varargin is not specified, return number of clusters
    if size(varargin, 2) == 0
        clusterData = numClusters;
        return;
    end;
    % Otherwise return points in cluster specified in varargin
    clusterId = varargin{1};
    idxSum = cumsum(idx);
    if clusterId == 1
        clustStart = 1;
    else
        clustStart = idxSum(clusterId - 1) + 1;
    end;
    clustEnd = idxSum(clusterId);
    clusterData = data(clustStart:clustEnd, :);
end;

%
% Helper function based on Ellipse_plot.m from Nima Moshtagh
% (nima@seas.upenn.edu)
% Copyright (c) 2009, Nima Moshtagh
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%
% Changes for plotClusters.m purpose:
% - Ellipses are plotted in black
% - Ellipse center is not plotted
% - Figure grid grid stays on.
%
function Ellipse_plot(A, C)

N = 20;


% check the dimension of the inputs: 2D or 3D
%--------------------------------------------
if length(C) == 3,
    Type = '3D';
elseif length(C) == 2,
    Type = '2D';
else
    display('Cannot plot an ellipse with more than 3 dimensions!!');
    return
end

% "singular value decomposition" to extract the orientation and the
% axes of the ellipsoid
[U D V] = svd(A);

if strcmp(Type, '2D'),
    % get the major and minor axes
    %------------------------------------
    a = 1/sqrt(D(1,1));
    b = 1/sqrt(D(2,2));

    theta = [0:1/N:2*pi+1/N];

    % Parametric equation of the ellipse
    %----------------------------------------
    state(1,:) = a*cos(theta); 
    state(2,:) = b*sin(theta);

    % Coordinate transform 
    %----------------------------------------
    X = V * state;
    X(1,:) = X(1,:) + C(1);
    X(2,:) = X(2,:) + C(2);
    
elseif strcmp(Type,'3D'),
    % generate the ellipsoid at (0,0,0)
    %----------------------------------
    a = 1/sqrt(D(1,1));
    b = 1/sqrt(D(2,2));
    c = 1/sqrt(D(3,3));
    [X,Y,Z] = ellipsoid(0,0,0,a,b,c,N);
    
    %  rotate and center the ellipsoid to the actual center point
    %------------------------------------------------------------
    XX = zeros(N+1,N+1);
    YY = zeros(N+1,N+1);
    ZZ = zeros(N+1,N+1);
    for k = 1:length(X),
        for j = 1:length(X),
            point = [X(k,j) Y(k,j) Z(k,j)]';
            P = V * point;
            XX(k,j) = P(1)+C(1);
            YY(k,j) = P(2)+C(2);
            ZZ(k,j) = P(3)+C(3);
        end
    end
end


% Plot the ellipse
%----------------------------------------
if strcmp(Type,'2D'),
    plot(X(1,:),X(2,:),'k');
else
    mesh(XX,YY,ZZ);
    axis equal
    hidden off
end

