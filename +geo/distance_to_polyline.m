function [d, v, it, p] = distance_to_polyline(line_pts, p0)
% Find the distance between given point(s) to a polyline
%
% INPUT
%   line_pts:       n*d, polyline points
%   p0:             m*d, the query points
%
% OUTPUT
%   d:              m*1, distance
%   v:              m*d, vector from query point to the nerest point on polyline
%   it:             m*2, line segment index and normalized length parameter t
%   v:              m*d, projected point on lines

m = size(p0, 1);
dim = size(p0, 2);

p = inputParser;
p.addRequired('pts', @(x) validateattributes(x, {'numeric'}, {'2d', 'ncols', dim}));
p.addRequired('p', @(x) validateattributes(x, {'numeric'}, {'2d'}));
p.parse(line_pts, p0);

d = nan(m, 1);
v = nan(m, dim);
it = nan(m, 2);
p = nan(m, dim);
for i = 1:m
    [tmp_dist, tmp_vn, tmp_t, tmp_p0] = p2line_dist(p0(i, :), line_pts(1:end - 1, :), line_pts(2:end, :));
    [d(i), tmp_idx] = min(tmp_dist);
    v(i, :) = tmp_vn(tmp_idx, :);
    it(i, :) = [tmp_idx, tmp_t(tmp_idx)];
    p(i, :) = tmp_p0(tmp_idx, :);
end
end

function [d, v, t, p0] = p2line_dist(p, xy1, xy2)
% Helper function. Find the distance between given point to several line segments
%
% INPUT
%   p:          1*d, query point
%   xy1, xy2:   n*d, several line segments
%
% OUTPUT
%   d:          scalar, distance
%   v:          n*d, vector from query point to every line segment
%   t:          n*1, normalized length parameter
%   p0:         n*d, projected point on lines

v1 = bsxfun(@minus, p, xy1);
v2 = xy2 - xy1;

t = min(max(sum(v1 .* v2, 2) ./ sum(v2.^2, 2), 0), 1);
p0 = bsxfun(@times, t, v2);
v = bsxfun(@minus, p0, v1);
d = sqrt(sum(v.^2, 2));
end
