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

% Let [u, s, v] = svd(jac), i.e. u * s * v' = jac.
% If rank(jac) = k, then there are only k values greater than 0 in s.
% We devide u, s, v into blocks,
% u = [u_o, u_n],     v = [v_o, v_n]
%      m*k  m*(m-k)        n*k  n*(n-k)
% Then jac = u_o * s(1:k, 1:k) * v_o'
% We will find that v_n is null space of jac: jac * v_n = u_o * s(1:k, 1:k) * v_o' * v_n = 0,
% and thus v_n is the tangent direction of contour.
% Let dx be in the orthogonal complement to v_n, i.e. perpendicular to contour. It can be
% expressed as dx = v_o * alpha. So,
%   |dy|^2 = dx' * jac' * jac * dx
%          = alpha' * v_o' * v_o * s_o * u_o' * u_o * s_o * v_o' * v_o * alpha
%          = alpha' * s_o^2 * alpha

h1 = max(std(data)) * 0.1;
rank_tol = 1e-10;

cq = cosd(0:5:360);
sq = sind(0:5:360);

hold on;
for i = 1:num
    [~, s, v] = svd(jac(1:3, :, i));
    jac_rank = sum(diag(s) > rank_tol);
    if jac_rank ~= 2
        continue;
    end
    s = max(s, s(1,1) * .01);
    dx = v(:, 1:2) * [cq / s(1, 1); sq / s(2, 2)] * h1;
    
    % Tolerance circle
    tmp_pts = bsxfun(@plus, data(i, :), dx');
    plot3(tmp_pts(:, sub_idx(1)), tmp_pts(:, sub_idx(2)), tmp_pts(:, sub_idx(3)), 'm', ...
        'linewidth', 1.5);
end
end