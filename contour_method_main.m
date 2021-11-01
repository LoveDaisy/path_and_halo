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

[r0, r0_ll, dr] = generate_healpix_grids(6);

r1 = nan(size(r0));
r2 = nan(size(r0));

entry_face_idx = 3;
exit_face_idx = 5;
entry_norm = face_norm(entry_face_idx, :);
exit_norm = face_norm(exit_face_idx, :);

valid_idx = r0 * entry_norm' < 0 & r0(:, 3) < 0;
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

seeds_idx = abs(target_bending - bending_angle) < dr;
checked_idx = false(size(seeds_idx));
[~, min_idx] = min(abs(target_bending - bending_angle));

[x, a, g_a] = find_bending_angle_solution(r0_ll(min_idx, :), target_bending, ...
    face_norm([entry_face_idx, exit_face_idx], :), [n, 1]);

h = 2;
max_h = 5;
min_h = 0.1;
h_factor = 0.05;
h_factor_min = 2e-4;
h_factor_max = 0.5;
res_num = 500;

res_store = nan(res_num, 4);  % [x1, x2, y]
res_store(1, :) = [x, a, h];
tan_g_store = nan(res_num, 2);
g_a = g_a * [0, -1; 1, 0];
tan_g_store(1, :) = g_a / norm(g_a);

i = 2;
% forward direction
non_shrink_cnt = 0;
while i <= res_num
    x0 = res_store(i-1, 1:2) + h * tan_g_store(i-1, :);
    [x, a, g_a] = find_bending_angle_solution(x0, target_bending, ...
        face_norm([entry_face_idx, exit_face_idx], :), [n, 1]);
    
    g_a = g_a * [0, -1; 1, 0];
    g_a = g_a / norm(g_a);
    
    % a simple adaptive schedule
    rho = norm(x - res_store(i-1, 1:2)) / abs(acos(dot(tan_g_store(i-1, :), g_a)));
    h = min(max(rho * h_factor, min_h), max_h);
    if isnan(a)
        non_shrink_cnt = 0;
        h_factor = h_factor / 2;
        h = res_store(i-1, 4) / 2;
        if h > min_h && h_factor > h_factor_min
            continue;
        else
            break;
        end
    elseif h < max_h && h_factor < h_factor_max
        non_shrink_cnt = non_shrink_cnt + 1;
        if  non_shrink_cnt > 5
            h_factor = h_factor * 1.2;
        end
    end
    res_store(i, :) = [x, a, h];
    tan_g_store(i, :) = g_a;
    
    % findout all checked seeds
    vec_a = bsxfun(@minus, x, r0_ll);
    vec_b = x - res_store(i-1, 1:2);
    t = sum(bsxfun(@times, vec_a, vec_b), 2) / norm(vec_b)^2;
    checked_idx = checked_idx | sqrt(sum((vec_a - t * vec_b).^2, 2)) < dr & (t >= 0 & t <= 1);
    
    i = i + 1;
end

% backword direction
if i <= res_num
    res_store(i, :) = res_store(1, :);
    tan_g_store(i, :) = tan_g_store(1, :);
    i = i + 1;
end
h_factor = 0.05;
h = 2;
non_shrink_cnt = 0;
while i <= res_num
    x0 = res_store(i-1, 1:2) - h * tan_g_store(i-1, :);
    [x, a, g_a] = find_bending_angle_solution(x0, target_bending, ...
        face_norm([entry_face_idx, exit_face_idx], :), [n, 1]);
    
    g_a = g_a * [0, -1; 1, 0];
    g_a = g_a / norm(g_a);
    
    % a simple adaptive schedule
    rho = norm(x - res_store(i-1, 1:2)) / abs(acos(dot(tan_g_store(i-1, :), g_a)));
    h = min(max(rho * h_factor, min_h), max_h);
    if isnan(a)
        non_shrink_cnt = 0;
        h_factor = h_factor / 2;
        h = res_store(i-1, 4) / 2;
        if h > min_h && h_factor > h_factor_min
            continue;
        else
            break;
        end
    elseif h < max_h && h_factor < h_factor_max
        non_shrink_cnt = non_shrink_cnt + 1;
        if  non_shrink_cnt > 5
            h_factor = h_factor * 1.2;
        end
    end
    res_store(i, :) = [x, a, h];
    tan_g_store(i, :) = g_a;
    
    % findout all checked seeds
    vec_a = bsxfun(@minus, x, r0_ll);
    vec_b = x - res_store(i-1, 1:2);
    t = sum(bsxfun(@times, vec_a, vec_b), 2) / norm(vec_b)^2;
    checked_idx = checked_idx | sqrt(sum((vec_a - t * vec_b).^2, 2)) < dr & (t >= 0 & t <= 1);
    
    i = i + 1;
end

%%
figure(1); clf;
plot(res_store(:, 1), res_store(:, 2), 'o');
axis equal; axis tight;
set(gca, 'xlim', [-180, 180], 'ylim', [-90, 90]);

hold on;
plot(r0_ll(seeds_idx, 1), r0_ll(seeds_idx, 2), 'ms');
plot(r0_ll(checked_idx, 1), r0_ll(checked_idx, 2), 'rx');