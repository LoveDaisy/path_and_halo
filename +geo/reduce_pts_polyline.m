function [pts, status] = reduce_pts_polyline(polyline, pts, varargin)
% Remove points that locate very close to given polyline.
% Result pts is normal formed. For llr representation, it means all longitudes will be in
% [-180, 180] (the same as atan2), all latitudes will be in [-90, 90] and rolls will be in [0, 360]
%
% INPUT
%   polyline:       n*d
%   pts:            m*d, all points to be checked.
%
% OPTION
%   'jac_fun':      function handle, Jacobian function, computes the Jacobian to a target space,
%                   in which distance can be evaluated. It must be in form of jac = fun(x) for
%                   x is n*d and jac is d1*d*n, where d is dimesion of input space and d1 is
%                   dimesion of output space.
%   'eps':          tolarence to decide whether a point is close to polyline.
%
% OUTPUT
%   pts:            m1*d, reduced points.
%   status:         struct

status.fun_eval_cnt = 0;
if isempty(pts)
    return;
end

p = inputParser;
p.addParameter('eps', 1e-4, @(x) validateattributes(x, {'double'}, {'positive'}));
p.addParameter('jac_fun', [], @(x) validateattributes(x, {'function_handle'}, {'scalar'}));
p.parse(varargin{:});

num = size(pts, 1);
dim1 = size(pts, 2);
if dim1 == 3
    pts0 = geo.normalize_llr(pts);
    polyline = geo.normalize_llr(polyline);
else
    pts0 = pts;
end

[d, v, itp] = geo.distance_to_polyline(polyline, pts0);
if ~isempty(p.Results.jac_fun)
    [~, jac] = p.Results.jac_fun(itp(:, 3:end));
    status.fun_eval_cnt = status.fun_eval_cnt + num;
    for i = 1:num
        d(i) = norm(jac(:, :, i) * v(i, :)');
    end
end

if dim1 == 4
    [d2, v2, itp2] = geo.distance_to_polyline(-polyline, pts0);
    if ~isempty(p.Results.jac_fun)
        [~, jac] = p.Results.jac_fun(itp2(:, 3:end));
        status.fun_eval_cnt = status.fun_eval_cnt + num;
        for i = 1:num
            d2(i) = norm(jac(:, :, i) * v2(i, :)');
        end
    end
    d = min(d, d2);
end

dup_idx = d < p.Results.eps;
pts = pts(~dup_idx, :);
end
