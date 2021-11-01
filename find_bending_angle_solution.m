function [x, a, g_a] = find_bending_angle_solution(x0, target, face_norm, n, varargin)
% INPUT
%    target:     m-vector or scalar, the target bending angle
%    x:          m*2, initial point
%    face_norm:  k*3, face normal
%    n:          k-vector, refractive index
% PARAEMTER (optional)
%   'eps':       scalar, default 1e-8
%   'MaxIter':   scalar, default 10

p = inputParser;
p.addParameter('eps', 1e-8);
p.addParameter('MaxIter', 5);
p.parse(varargin{:});

x = x0;
max_step = 50;

[~, a, ~, g_a] = bending_angle_with_gradient(x, face_norm, n);
if isnan(a)
    x = nan(size(x0));
    g_a(:) = nan;
    return
end

da = target - a;
iter_num = 1;
while abs(da) > p.Results.eps && iter_num < p.Results.MaxIter
    step = da / norm(g_a)^2;
    step = min(abs(step), max_step) * sign(step);
    dx = step * g_a;
    
    alpha = 2;
    a = nan;
    while (isnan(a) || abs(a - target) > abs(da) * 0.8) && alpha > 0.1
        alpha = alpha / 2;
        [~, a, ~, g_a] = bending_angle_with_gradient(x + dx * alpha, face_norm, n);
    end
    
    x = x + dx * alpha;
    da = target - a;
    iter_num = iter_num + 1;
end

if abs(da) > p.Results.eps
    x(:) = nan;
    a = nan;
    g_a(:) = nan;
end
end
