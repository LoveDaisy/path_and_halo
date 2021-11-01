clear; close all; clc;

face_norm = [0, 0, 1;     % face 1
    0, 0, -1;             % face 2
    -sqrt(3)/2, 1/2, 0;   % face 3
    0, 1, 0;              % face 4
    sqrt(3)/2, 1/2, 0;    % face 5
    sqrt(3)/2, -1/2, 0;   % face 6
    0, -1, 0;             % face 7
    -sqrt(3)/2, -1/2, 0]; % face 8
n = 1.31;

[r0, r0_ll, dr] = generate_healpix_grids(5);

r1 = nan(size(r0));
r2 = nan(size(r0));

entry_face_idx = 3;
exit_face_idx = 5;
entry_norm = face_norm(entry_face_idx, :);
exit_norm = face_norm(exit_face_idx, :);

valid_idx = r0 * entry_norm' < 0;
r1(valid_idx, :) = refract_with_gradient(r0(valid_idx, :), entry_norm, 1, n);
valid_idx = valid_idx & (r1 * exit_norm' > 0);
r2(valid_idx, :) = refract_with_gradient(r1(valid_idx, :), exit_norm, n, 1);
valid_idx = valid_idx & sum(r2.^2, 2) > 1e-4;
r2(~valid_idx, :) = nan;

bending_angle = acosd(sum(r0 .* r2, 2));
bending_angle_max = max(bending_angle);
bending_angle_min = min(bending_angle);

%%
target_bending = 43;

seeds_idx = abs(target_bending - bending_angle) < dr * 2;
checked_idx = false(size(seeds_idx));
[~, min_idx] = min(abs(target_bending - bending_angle));

[x, a, g_a] = find_bending_angle_solution(r0_ll(min_idx, :), target_bending, ...
    face_norm([entry_face_idx, exit_face_idx], :), [n, 1]);

h = 2;
max_h = 5;
min_h = 0.001;
h_factor = 0.05;
h_factor_min = 1e-4;
h_factor_max = 0.5;
hh = -0.1;

res_num = 100;
res_store = nan(res_num, 5);  % [x1, x2, y]
res_store(1, :) = [x, a, h, 0];
grad_store = nan(res_num, 4);  % [original_gradient, normalized_co_gradient]

grad_store(1, 1:2) = g_a;
g_a = g_a * [0, -1; 1, 0];
grad_store(1, 3:4) = g_a / norm(g_a);

i = 2;
% forward direction
non_shrink_cnt = 0;
while i <= res_num
    if i <= 3
        x0 = res_store(i-1, 1:2) + h * grad_store(i-1, 3:4);
    else
        % use 2nd order prediction
        tmp_t = cumsum(res_store(i-3:i-1, 5));
        tmp_x = res_store(i-3:i-1, 1:2);
        new_t = tmp_t(3) + h;
        new_x = (new_t - tmp_t(1)) * (new_t - tmp_t(2)) * (tmp_t(1) - tmp_t(2)) * tmp_x(3, :) + ...
            (new_t - tmp_t(2)) * (new_t - tmp_t(3)) * (tmp_t(2) - tmp_t(3)) * tmp_x(1, :) + ...
            (new_t - tmp_t(3)) * (new_t - tmp_t(1)) * (tmp_t(3) - tmp_t(1)) * tmp_x(2, :);
        x0 = -new_x / ((tmp_t(1) - tmp_t(2)) * (tmp_t(2) - tmp_t(3)) * (tmp_t(3) - tmp_t(1)));
    end
    x0 = x0 + hh * dr * grad_store(i-1, 1:2) / norm(grad_store(i-1, 1:2));
    [x, a, g_a] = find_bending_angle_solution(x0, target_bending, ...
        face_norm([entry_face_idx, exit_face_idx], :), [n, 1]);
    if i >= 3 && dot(x - res_store(i-1, 1:2), grad_store(i-1, 3:4)) / ...
            norm(x - res_store(i-1, 1:2)) < 0
        break;
    end
    
    grad_store(i, 1:2) = g_a;
    g_a = g_a * [0, -1; 1, 0];
    g_a = g_a / norm(g_a);
    
    % a simple adaptive schedule
    if isnan(a)
        non_shrink_cnt = 0;
        h = h / 2;
        if h > min_h
            continue;
        else
            break;
        end
    elseif h < max_h
        non_shrink_cnt = non_shrink_cnt + 1;
        if  non_shrink_cnt > 3
            h = h * 1.2;
        end
    end
    s = norm(x - res_store(i-1, 1:2));
    res_store(i, :) = [x, a, h, s];
    grad_store(i, 3:4) = g_a;
    
    % findout all checked seeds
    vec_a = bsxfun(@minus, x, r0_ll);
    vec_b = x - res_store(i-1, 1:2);
    t = sum(bsxfun(@times, vec_a, vec_b), 2) / norm(vec_b)^2;
    checked_idx = checked_idx | sqrt(sum((vec_a - t * vec_b).^2, 2)) < dr & (t >= 0 & t <= 1);
    
    i = i + 1;
end

% backword direction
start_i = i;
if i <= res_num
    res_store(i, :) = res_store(1, :);
    grad_store(i, :) = grad_store(1, :);
    i = i + 1;
end
h = 2;
non_shrink_cnt = 0;
while i <= res_num
    if i <= start_i + 2
        x0 = res_store(i-1, 1:2) - h * grad_store(i-1, 3:4);
    else
        % use 2nd order prediction
        tmp_t = cumsum(-res_store(i-3:i-1, 5));
        tmp_x = res_store(i-3:i-1, 1:2);
        new_t = tmp_t(3) - h;
        new_x = (new_t - tmp_t(1)) * (new_t - tmp_t(2)) * (tmp_t(1) - tmp_t(2)) * tmp_x(3, :) + ...
            (new_t - tmp_t(2)) * (new_t - tmp_t(3)) * (tmp_t(2) - tmp_t(3)) * tmp_x(1, :) + ...
            (new_t - tmp_t(3)) * (new_t - tmp_t(1)) * (tmp_t(3) - tmp_t(1)) * tmp_x(2, :);
        x0 = -new_x / ((tmp_t(1) - tmp_t(2)) * (tmp_t(2) - tmp_t(3)) * (tmp_t(3) - tmp_t(1)));
    end
    x0 = x0 + hh * dr * grad_store(i-1, 1:2) / norm(grad_store(i-1, 1:2));
    [x, a, g_a] = find_bending_angle_solution(x0, target_bending, ...
        face_norm([entry_face_idx, exit_face_idx], :), [n, 1]);
    if i >= 3 && dot(x - res_store(i-1, 1:2), grad_store(i-1, 3:4)) / ...
            norm(x - res_store(i-1, 1:2)) > 0
        break;
    end
    
    grad_store(i, 1:2) = g_a;
    g_a = g_a * [0, -1; 1, 0];
    g_a = g_a / norm(g_a);
    
    % a simple adaptive schedule
    if isnan(a)
        non_shrink_cnt = 0;
        h = h / 2;
        if h > min_h
            continue;
        else
            break;
        end
    elseif h < max_h
        non_shrink_cnt = non_shrink_cnt + 1;
        if  non_shrink_cnt > 3
            h = h * 1.2;
        end
    end
    s = norm(x - res_store(i-1, 1:2));
    res_store(i, :) = [x, a, h, s];
    grad_store(i, 3:4) = g_a;
    
    % findout all checked seeds
    vec_a = bsxfun(@minus, x, r0_ll);
    vec_b = x - res_store(i-1, 1:2);
    t = sum(bsxfun(@times, vec_a, vec_b), 2) / norm(vec_b)^2;
    checked_idx = checked_idx | sqrt(sum((vec_a - t * vec_b).^2, 2)) < dr & (t >= 0 & t <= 1);
    
    i = i + 1;
end

%%
figure(1); clf;
plot(res_store(:, 1), res_store(:, 2), '-o');
axis equal; axis tight;
set(gca, 'xlim', [-180, 180], 'ylim', [-90, 90]);

hold on;
plot(r0_ll(seeds_idx & bending_angle < target_bending, 1), ...
    r0_ll(seeds_idx & bending_angle < target_bending, 2), 'ks');
plot(r0_ll(seeds_idx & bending_angle > target_bending, 1), ...
    r0_ll(seeds_idx & bending_angle > target_bending, 2), 'ms');
plot(r0_ll(checked_idx, 1), r0_ll(checked_idx, 2), 'rx');