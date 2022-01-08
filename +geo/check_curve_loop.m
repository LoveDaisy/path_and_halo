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
%   'eps':      scalar, point-to-line distance tolerance
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
p.parse(x, varargin{:});

[d, ~, ~] = geo.distance_to_polyline(x(1:end - 1, :), x(end, :));
last_x = [];
while d < p.Results.eps
    last_x = x(end, :);
    x = x(1:end - 1, :);
    [d, ~, ~] = geo.distance_to_polyline(x(1:end - 1, :), x(end, :));
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
    v1 = x(1, :) - last_x;
    v2 = x(2, :) - x(1, :);
    closed = dot(v1, v2) / norm(v1) / norm(v2) > 0.3;
end

end
