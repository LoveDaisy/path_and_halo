function [area_factor, vtx] = project_polygon_intersection(vtx_from, vtx_to, ray)
% Project a polygon onto another one and find out the intersection polygon and the area factor.
% The area factor is defined as a ratio of the area of intersection polygon to the area of original
% polygon.
%
% INPUT
%   vtx_from:   n*d, vertex of original polygon. d > 2
%   vtx_to:     m*d, vertex of target polygon. it must be *CONVEX*
%   ray:        1*d, projection ray
%
% OUTPUT
%   area_factor:    scalar, the ratio that projected intersection area to original area.
%   vtx :           k*d, vertex of intersetion polygon. Clearly it is coplanar to vtx_to.

area_factor = 0;
vtx = vtx_from;

% Check dimesion
dim = size(vtx_from, 2);
if dim <= 2 || size(vtx_to, 2) ~= dim || size(ray, 2) ~= dim
    error('dimension less than 2 or dimension mismatch!');
end

if size(vtx_to, 1) < 3
    error('vtx_to must be more than three points!');
end

% Check if all vertexes are coplanar
if ~geo.check_coplanar(vtx_from) || ~geo.check_coplanar(vtx_to)
    warning('vtx_from and vtx_to must be coplanar!');
    return;
end

% Say original point is p0, projection ray is r, target polygon is qi, i = 1, 2, ...
% and projection of p0 is pq, then we will see:
% (pq - q1) - (p0 - q1) = t * r  ......(1)
% In order to compute area we need to express pq in 2D form:
% pq - q1 = [q2 - q1, q3 - q1] * uv = Bq * uv, where Bq is d*2 and uv is 2*1. So (1) becomes
% Bq * uv - p0 + q1 = t * r, or [Bq, -r] * [uv; t] = p0 - q1  ......(2)
Bq = [vtx_to(2, :) - vtx_to(1, :); vtx_to(3, :) - vtx_to(1, :)];
uvt = bsxfun(@minus, vtx_from, vtx_to(1, :)) / [Bq; -ray];
uv_from = uvt(:, 1:2);
uv_to = bsxfun(@minus, vtx_to, vtx_to(1, :)) / Bq;

% Find the intersection.
uv_int = geo.polygon2d_intersection(uv_from, uv_to);
vtx = bsxfun(@plus, uv_int * Bq, vtx_to(1, :));
if size(vtx, 1) < 3
    return;
end

area_from = sum(uv_from(:, 1) .* [uv_from(2:end, 2); uv_from(1, 2)] - ...
    uv_from(:, 2) .* [uv_from(2:end, 1); uv_from(1, 1)]);
area_int = sum(uv_int(:, 1) .* [uv_int(2:end, 2); uv_int(1, 2)] - ...
    uv_int(:, 2) .* [uv_int(2:end, 1); uv_int(1, 1)]);
area_factor = area_int / area_from;
end
