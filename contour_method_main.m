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
sun_altitude = 20;
sun_longitude = 180;
ray_in = -[cosd(sun_altitude) * cosd(sun_longitude), cosd(sun_altitude) * sind(sun_longitude), ...
    sind(sun_altitude)];

target_bending = 21.8395;

[~, min_idx] = min(abs(target_bending - bending_angle));
start_pt1 = r0_ll(min_idx, :);
start_pt2 = start_pt1; start_pt2(2) = -start_pt1(2);
[x, a, g_a] = find_bending_angle_solution(start_pt1, target_bending, ...
    face_norm([entry_face_idx, exit_face_idx], :), [n, 1]);

res_num = 200;
res_store = nan(res_num, 3);  % [x1, x2, y]
res_store(1, :) = [x, a];
grad_store = nan(res_num, 2);
grad_store(1, :) = g_a;

h = 5;
i = 2;
k = 1;
config = [start_pt1, 1; start_pt1, -1; start_pt2, 1; start_pt2, -1];
while i <= res_num
    x0 = res_store(i-1, 1:2) + h * [g_a(2), -g_a(1)] * config(k, 3);
    [x, a, g_a] = find_bending_angle_solution(x0, target_bending, ...
        face_norm([entry_face_idx, exit_face_idx], :), [n, 1]);
    change_config = isnan(a);
    if change_config && k < 4
        k = k + 1;
        [x, a, g_a] = find_bending_angle_solution(config(k, 1:2), target_bending, ...
            face_norm([entry_face_idx, exit_face_idx], :), [n, 1]);
    end
    if isnan(a)
        break;
    end
    res_store(i, :) = [x, a];
    grad_store(i, :) = g_a;
    
    % a simple adaptive schedule
    dv_norm = acosd(dot(grad_store(i-1, :), g_a) / norm(g_a) / norm(grad_store(i-1, :)));
    if ~change_config && dv_norm < 1.5 && h < 10
        h = h * 1.5;
    elseif ~change_config && dv_norm > 3 && h > 0.01
        h = h / 1.5;
    else
        i = i + 1;
    end
end

figure(1); clf;
plot(res_store(:, 1), res_store(:, 2), 'o');
axis equal; axis tight;
set(gca, 'xlim', [-180, 180], 'ylim', [-90, 90]);