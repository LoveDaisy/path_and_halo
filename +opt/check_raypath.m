function res = check_raypath(crystal, raypath)
% Check if a raypath is valid for a given crystal.
%
% INPUT
%   crystal:        a struct
%   raypath:        n*1, face number indicating a raypath
%
% OUTPUT
%   res:            true or false

num = length(raypath);
res = true;
if num == 1
    return;
end

% A point in a polygon can be expressed as: p = sum(v_i * x_i), s.t. sum(x_i) = 1 and x_i >= 0;
% So if there is a ray go through all polygons, then it must be hold for all k = 2, 3, ..., N-1
%   P_1 * (s_k * alpha) + P_N * (t_k * beta) = P_k * gamma_k
%   ^~~~         ^~~~~    ^~~~         ^~~~    ^~~~  ^~~~~~
%    d*n1         n1*1     d*nN         nN*1    d*nk  nk*1
% where P_k is n_k vertexes for polygon k.
%
% We can thus construct an optimization problem:
%   min  sum_k |P_1 * (s_k * alpha) + P_N * (t_k * beta) - P_k * gamma_k|^2
%   s.t. 0 <= alpha, beta, gamma_k, s_k, t_k <= 1
%        sum(alpha) = 1, sum(beta) = 1, sum(gamma_k) = 1, s_k + t_k = 1
% If the minmumn is 0 then the solution exists.
%
% This solution doesn't count in refraction. In a corner case, a ray goes through all polygons,
% but its incident angle is too large that it cannot go out of the last surface (total reflection
% occurs). So we have to find as small incident angle as possible, and check if there is total
% reflection.
%
% We then modify the optimization problem:
%   min  max(-dot(normal_N, r), dot(normal_1, r)) + sum_k( |r_k|^2 ) * factor
%   s.t. r_k = P_1 * (s_k * alpha) + P_N * (t_k * beta) - P_k * gamma_k
%        r = normalize(-P_1 * (s_k * alpha) + P_N * (t_k * beta))
%        0 <= alpha, beta, gamma_k, s_k, t_k <= 1
%        sum(alpha) = 1, sum(beta) = 1, sum(gamma_k) = 1, s_k + t_k = 1

entry_exit_norm = nan(2, 3);
entry_exit_norm(1, :) = crystal.face_norm(raypath(1), :);
polygon_list = cell(num, 1);
polygon_list{1} = crystal.vtx(crystal.face{raypath(1)}, :);

curr_vtx = crystal.vtx;
curr_norm = crystal.face_norm;
for i = 2:num-1
    p = curr_vtx(crystal.face{raypath(i)}, :);
    p0 = mean(p);
    m = (eye(3) - 2 * curr_norm(raypath(i), :)' * curr_norm(raypath(i), :));
    next_vtx = bsxfun(@plus, bsxfun(@minus, curr_vtx, p0) * m, p0);
    next_norm = curr_norm * m;

    polygon_list{i} = p;
    curr_vtx = next_vtx;
    curr_norm = next_norm;
end

polygon_list{num} = curr_vtx(crystal.face{raypath(end)}, :);
entry_exit_norm(2, :) = curr_norm(raypath(end), :);

% Find each polygon center as start value for optimization
polygon_centers = zeros(num, 3);
for i = 1:num
    polygon_centers(i, :) = mean(polygon_list{i});
end
d = sqrt(sum(diff(polygon_centers).^2, 2));
t = [0; cumsum(d)] / sum(d);
s = 1 - t;

nk = zeros(num, 1);
for i = 1:num
    nk(i) = size(polygon_list{i}, 1);
end

% Set initial value
x_alpha = ones(nk(1), 1) / nk(1);
x_beta = ones(nk(end), 1) / nk(end);
x_gamma = zeros(sum(nk(2:num-1)), 1);
offset = 0;
for i = 2:num-1
    x_gamma(offset+1:offset+nk(i)) = ones(nk(i), 1) / nk(i);
    offset = offset + nk(i);
end
x0 = [x_alpha; x_beta; x_gamma; s(2:num-1); t(2:num-1)];

% Make Aeq for optimization
Aeq = zeros(num, length(x0));
Aeq(1, 1:nk(1)) = 1;
Aeq(2, (1:nk(end))+nk(1)) = 1;
offset = nk(1) + nk(end);
for i = 2:num-1
    Aeq(i+1, offset+(1:nk(i))) = 1;
    offset = offset + nk(i);
end
for i = 2:num-1
    Aeq(num-1+i, offset+[1, num-1]) = 1;
    offset = offset + 1;
end
beq = ones(size(Aeq, 1), 1);

edge_eps = 1e-4;
lb = zeros(length(x0), 1) + edge_eps;
ub = ones(length(x0), 1) - edge_eps;

res = false;

% Stage 1: ignore incident angle at entry and exit surface
option = optimoptions('fmincon', 'Display', 'off', 'MaxFunctionEvaluations', 3000, ...
    'ObjectiveLimit', 0.5);

obj_fun = @(x) solve_prog(x, polygon_list, [], edge_eps * 1.05);
val0 = obj_fun(x0);
if isnan(val0)
    return;
end

[x1, ~, opt_flag, opt_output] = fmincon(obj_fun, x0, [], [], Aeq, beq, lb, ub, ...
    [], option);

[val1, ~, rk_loss, ~, is_edge] = obj_fun(x1);
cos_angle = find_incident_angle(x1, polygon_list, entry_exit_norm);
res = (opt_flag >= 0 || opt_flag == -3) && all(~is_edge) && rk_loss < 1e-6;
if (res && cos_angle > cos(asin(1 / crystal.n))) || ~res
    return;
end

% Stage 2: find out minimun incident angle
option = optimoptions(option, 'ObjectiveLimit', -0.7);
obj_fun = @(x) solve_prog(x, polygon_list, entry_exit_norm, edge_eps * 1.05);

[x2, ~, opt_flag, opt_output] = fmincon(obj_fun, x1, [], [], Aeq, beq, lb, ub, ...
    [], option);

[val2, neg_cos, rk_loss, ~, is_edge] = obj_fun(x2);
res = (opt_flag >= 0 || opt_flag == -3) && all(~is_edge) && ...
    -neg_cos > cos(asin(1 / crystal.n)) && rk_loss < 1e-6;
end

function [val, neg_cos, rk_loss, x_loss, is_edge] = solve_prog(x, polys, entry_exit_norm, edge_eps)
% x:  [alpha, beta, gamma_2, gamma_3, ..., gamma_N-1, s_2, s_3, ..., s_N-1]
%      n1     nN     n2       n3             nN-1      1    1          1
% The optimization problem is:
%   min  max(-dot(normal_N, r), dot(normal_1, r)) + sum_k( |r_k|^2 ) * factor
%   s.t. r_k = P_1 * (s_k * alpha) + P_N * ((1 - s_k) * beta) - P_k * gamma_k
%        r = normalize(-P_1 * alpha + P_N * beta)
%        0 <= alpha, beta, gamma_k, s_k <= 1
%        sum(alpha) = 1, sum(beta) = 1, sum(gamma_k) = 1

num = length(polys);
nk = zeros(num, 1);
for i = 1:num
    nk(i) = size(polys{i}, 1);
end
total_k_n = sum(nk(2:num-1));

offset = 0;
alpha_ = x((1:nk(1))+offset);
offset = offset + nk(1);
beta_ = x((1:nk(end))+offset);
offset = offset + nk(end);
gamma_ = x((1:total_k_n)+offset);
offset = offset + total_k_n;
s_ = x((1:num-2)+offset);
offset = offset + num - 2;
t_ = x((1:num-2)+offset);

is_edge = false(num, 1);
is_edge(1) = sum(alpha_ >= edge_eps) < 3;
is_edge(2) = sum(beta_ >= edge_eps) < 3;

offset = 0;
rk = zeros(num-2, 3);
for i = 2:num-1
    rk(i-1, :) = s_(i-1) * alpha_' * polys{1} + t_(i-1) * beta_' * polys{end} - ...
        gamma_(offset+1:offset+nk(i))' * polys{i};
    is_edge(i+1) = sum(gamma_(offset+1:offset+nk(i)) >= edge_eps) < 3;
    offset = offset + nk(i);
end

x_loss = sum([alpha_; beta_; gamma_].^2);
if num > 2
    rk_loss = mean(rk(:).^2);
else
    rk_loss = 0;
end

if ~isempty(entry_exit_norm)
    neg_cos = -find_incident_angle(x, polys, entry_exit_norm);
else
    neg_cos = 0;
end
val = neg_cos + rk_loss * 1e8 + x_loss * 1e-3;
end

function cos_angle = find_incident_angle(x, polys, entry_exit_norm)
num = length(polys);
nk = zeros(num, 1);
for i = 1:num
    nk(i) = size(polys{i}, 1);
end

offset = 0;
alpha_ = x((1:nk(1))+offset);
offset = offset + nk(1);
beta_ = x((1:nk(end))+offset);

r = -alpha_' * polys{1} + beta_' * polys{end};
r = r / norm(r);

cos_angle = min([-dot(entry_exit_norm(1, :), r), dot(entry_exit_norm(2, :), r)]);
end
