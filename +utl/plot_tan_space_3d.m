function plot_tan_space_3d(data, fdf, sub_idx, varargin)
% Plot tangent vector & (unit) tolerance circle for data points.
% Let J be Jacobian of function fdf, which is an m*n matrix. Thus dy = J * dx.
% m rows (assume they are linearly independent) of J span an m-subspace in R^n, and its
% null space is R^(n - m).
% The direction of contour (keep a constant output value y) lies in null space of J, and the
% tolerance circle lies in row space of J.
% Let dx be a linear combination of rows of J, dx = J' * a, so dy = J * dx = J * J' * a
% Let M = J * J' and clearly M is symmetrical. Then dy = M * a, and |dy|^2 = a' * M' * M * a
% For tolerance circle we have |dy|^2 == const. We make orthogonal diagonalization of M as
% M = Q' * S * Q, and we get |dy|^2 = a' * Q' * S^2 * Q * a
% Let b = Q * a, then |dy|^2 = b' * diag([s1^2, s2^2]) * b
% If b = [1/s1 cos(q); 1/s2 sin(q)] or b = [1/s1 sin(q); 1/s2 cos(q)], we get |dy|^2 == 1
%
% INPUT
%   data:       n*d
%   fdf:        function handle
%   sub_idx:    length-3-vector, indicating which three sub-idx of data will be displayed

if isempty(sub_idx)
    sub_idx = [1, 2, 3];
end

num = size(data, 1);
[~, jac] = fdf(data);

h1 = max(std(data)) * 5;
h2 = -inf;
for i = 1:num
    h2 = max(h2, norm(jac(1, :, i)));
    h2 = max(h2, norm(jac(2, :, i)));
end
h2 = 0.04 / h2 * h1;

cq = cosd(0:5:360);
sq = sind(0:5:360);

hold on;
for i = 1:num
    jjt = jac(1:2, :, i) * jac(1:2, :, i)';
    [u, s, v] = svd(jjt);
    s = s + eye(2) * max(s(:)) * 1e-4;
    if det(u) > 0
        matQt = u;
    elseif det(v) > 0
        matQt = v;
    else
        matQt = u * diag([1, -1]);
    end
    a = matQt * [1/s(1, 1) * cq; 1/s(2, 2) * sq];
    dx = a' * jac(1:2, :, i) * h1;
    
    % Tolerance circle
    tmp_pts = bsxfun(@plus, data(i, :), dx);
    plot3(tmp_pts(:, sub_idx(1)), tmp_pts(:, sub_idx(2)), tmp_pts(:, sub_idx(3)), 'm', ...
        'linewidth', 1.5);
    
    % Tangent vector
    plot3([0, jac(1, sub_idx(1), i)] * h2 + data(i, sub_idx(1)), ...
        [0, jac(1, sub_idx(2), i)] * h2 + data(i, sub_idx(2)), ...
        [0, jac(1, sub_idx(3), i)] * h2 + data(i, sub_idx(3)), 'k', 'linewidth', 1.2);
    plot3([0, jac(2, sub_idx(1), i)] * h2 + data(i, sub_idx(1)), ...
        [0, jac(2, sub_idx(2), i)] * h2 + data(i, sub_idx(2)), ...
        [0, jac(2, sub_idx(3), i)] * h2 + data(i, sub_idx(3)), 'k', 'linewidth', 1.2);
end
end