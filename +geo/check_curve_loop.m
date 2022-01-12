function [closed, x] = check_curve_loop(x, varargin)
% Check whether a curve is a closed loop. And reduce redundant points (from the tail) if neccessary.
% If one end lies in the middle of a line segment of the curve, it is closed;
% If distance from one end to the other end is less than either segment length, it is closed.
% It is NOT closed if there is a self-cross.
% Assume there is no sharp angle on curve.
%
% INPUT
%   x:          n*d, curve points
%
% OPTION
%   'eps':      scalar, point-to-line distance tolerance. Default is 1e-3.
%   'IntStep':  scalar, step used in interpolation. If greater than 0, input polyline will be
%               interpolated by this step. Ignore when it is less than 0. Default is -1;
%
% OUTPUT
%   closed:     bool
%   x:          m*d, reduced curve points

num = size(x, 1);
if num < 3
    closed = false;
    return;
end

p = inputParser;
p.addRequired('x', @(x) validateattributes(x, {'numeric'}, {'2d'}));
p.addParameter('eps', 1e-3, @(x) validateattributes(x, {'double'}, {'positive'}));
p.addParameter('IntStep', -1, @(x) validateattributes(x, {'double'}, {'real'}));
p.parse(x, varargin{:});

if p.Results.IntStep > 0 && size(x, 1) > 2
    polyline = geo.interp_curve(x(1:end - 1, :), p.Results.IntStep);
else
    polyline = x(1:end - 1, :);
end
[d, ~, ~] = geo.distance_to_polyline(polyline, x(end, :));
last_x = [];
while d < p.Results.eps
    last_x = x(end, :);
    x = x(1:end - 1, :);
    if p.Results.IntStep > 0 && size(x, 1) > 2
        polyline = geo.interp_curve(x(1:end - 1, :), p.Results.IntStep);
    else
        polyline = x(1:end - 1, :);
    end
    [d, ~, ~] = geo.distance_to_polyline(polyline, x(end, :));
end

num = size(x, 1);
if num < 3
    closed = false;
    return;
end

if isempty(last_x)
    d0 = norm(x(end, :) - x(1, :));
    d1 = norm(x(1, :) - x(2, :));
    dn = norm(x(end - 1, :) - x(end, :));
    closed = d0 < d1 || d0 < dn;
else
    v1 = x(1, :) - x(end, :);
    v2 = x(2, :) - x(1, :);
    closed = dot(v1, v2) / norm(v1) / norm(v2) > 0.3;
end

end
