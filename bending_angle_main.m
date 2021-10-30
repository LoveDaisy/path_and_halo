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

n_side = 2^5;
n_pix = nSide2nPix(n_side);
dr = sqrt(4 * pi / n_pix) * 180 / pi * 0.5;

[ray_in_x, ray_in_y, ray_in_z] = pix2vec(n_side, 1:n_pix);
r0 = -[ray_in_x', ray_in_y', ray_in_z'];
clear ray_in_x ray_in_y ray_in_z;
r0_ll = [atan2d(r0(:, 2), r0(:, 1)), asind(r0(:, 3) ./ sqrt(sum(r0.^2, 2)))];  % longitude, latitude

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

% lon = -32;
% lat = -24;
% ray_out = [cosd(lat) * cosd(lon), cosd(lat) * sind(lon), sind(lat)];
% target_bending = acosd(dot(ray_out, ray_in));
target_bending = 38;

[~, min_idx] = min(abs(target_bending - bending_angle));
x = r0_ll(min_idx, :);
[~, a, ~, g_angle] = ...
    bending_angle_with_gradient(x, face_norm([entry_face_idx, exit_face_idx], :), [n, 1]);
da = target_bending - a;
while abs(da) > 1e-8
    x = da / norm(g_angle)^2 * g_angle + x;
    [~, a, ~, g_angle] = ...
        bending_angle_with_gradient(x, face_norm([entry_face_idx, exit_face_idx], :), [n, 1]);
    da = target_bending - a;
end

res_num = 500;
res_store = zeros(res_num + 1, 3);  % [x1, x2, y]
res_store(1, :) = [x, a];
h = 1;
for i = 1:res_num
    x = res_store(i, 1:2);
    [~, a, ~, g_angle] = ...
        bending_angle_with_gradient(x, face_norm([entry_face_idx, exit_face_idx], :), [n, 1]);
    da = target_bending - a;
    x = da / norm(g_angle)^2 * g_angle + x;
    [~, a, ~, g_angle] = ...
        bending_angle_with_gradient(x, face_norm([entry_face_idx, exit_face_idx], :), [n, 1]);
    res_store(i, :) = [x, a];
    t_a = [g_angle(2), -g_angle(1)];
    res_store(i+1, 1:2) = x + t_a * h;
end

figure(1); clf;
plot(res_store(:, 1), res_store(:, 2), 'o');
axis equal; axis tight;
set(gca, 'xlim', [-180, 180], 'ylim', [-90, 90]);