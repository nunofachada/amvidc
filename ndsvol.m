function v = ndsvol(n, r)
% NDSVOL Determine the volume of a n-dimensional sphere.
% 
%   v = NDSVOL(n, r)
%
% Parameters:
%    n - Number of dimensions.
%    r - Radius.
% Output:
%    v - Volume of a n-dimensional sphere with radius = r.
%

%  N. Fachada
%  Instituto Superior Técnico, Lisboa, Portugal

% Use known unit sphere volume to speedup computation
switch n
    case 2
        % 2D case (area of circle)
        v = pi * r^2;
    case 3
        % 3D case (volume of sphere)
        v = 4 / 3 * pi * r^3;
    otherwise
        % ND case, as described here:
        % http://mathworld.wolfram.com/Hypersphere.html
        v = (2 * pi^(n / 2) * r^n) / (n * gamma(n / 2));
end;
