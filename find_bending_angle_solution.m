function [x, a, g_a] = find_bending_angle_solution(x0, target, face_norm, n, varargin)
% INPUT
%    target:      m-vector or scalar, the target bending angle
%    x:          m*2, initial point
%    face_norm:  k*3, face normal
%    n:          k-vector, refractive index

p = inputParser;
p.addParameter('eps', 1e-8);
p.addParameter('MaxIter', 10);
p.parse(varargin{:});

x = x0;
[~, a, ~, g_a] = bending_angle_with_gradient(x, face_norm, n);
da = target - a;
iter_num = 1;
while abs(da) > p.Results.eps && iter_num < p.Results.MaxIter
    dx = da / norm(g_a)^2 * g_a;
    x = x + dx;
    [~, a, ~, g_a] = bending_angle_with_gradient(x, face_norm, n);
    da = target - a;
    iter_num = iter_num + 1;
end

if iter_num >= p.Results.MaxIter
    warning('Iteration exceed limit %d!\n', p.Results.MaxIter);
end
end
