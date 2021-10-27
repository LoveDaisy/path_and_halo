clear; close all; clc;

n = 1.31;
norm1 = [-sqrt(3)/2, 1/2, 0];  % face 3
norm2 = [sqrt(3)/2, 1/2, 0];   % face 5

n_side = 2^6;
n_pix = nSide2nPix(n_side);
dr = sqrt(4 * pi / n_pix) * 180 / pi * 0.5;

[ray_in_x, ray_in_y, ray_in_z] = pix2vec(n_side, 1:n_pix);
r0 = -[ray_in_x', ray_in_y', ray_in_z'];
clear ray_in_x ray_in_y ray_in_z;
r0_ll = [atan2d(r0(:, 2), r0(:, 1)), asind(r0(:, 3) ./ sqrt(sum(r0.^2, 2)))];  % longitude, latitude

r1 = nan(size(r0));
r2 = nan(size(r0));

valid_idx = r0 * norm1' < 0;
r1(valid_idx, :) = refract(r0(valid_idx, :), norm1, 1, n);
valid_idx = valid_idx & (r1 * norm2' > 0);
r2(valid_idx, :) = refract(r1(valid_idx, :), norm2, n, 1);
valid_idx = valid_idx & sum(r2.^2, 2) > 1e-4;
r2(~valid_idx, :) = nan;

bending_angle = acosd(sum(r0 .* r2, 2));
bending_angle_max = max(bending_angle);
bending_angle_min = min(bending_angle);

%%
% figure(1); clf;
% scatter(atan2d(r0(:, 2), r0(:, 1)), asind(r0(:, 3) ./ sqrt(sum(r0.^2, 2))), ...
%     10, bending_angle, 'filled');
% axis equal; axis tight;
% box on;

%%
sun_altitude = 10;  % degree
ray_in = -[cosd(sun_altitude) * cosd(-90), cosd(sun_altitude) * sind(-90), sind(sun_altitude)];
ray_out = -r0;

crystal_zenith_mean = 90;
crystal_zenith_std = 1;

halo_img = zeros(n_pix, 1);
for i = 1:n_pix
    curr_ray_out = ray_out(i, :);
    curr_bending = acosd(dot(curr_ray_out, ray_in));
    if curr_bending < bending_angle_min || curr_bending > bending_angle_max
        continue;
    end
    
    curr_idx = abs(bending_angle - curr_bending) < dr & valid_idx;
    if ~any(curr_idx)
        continue;
    end
    
    curr_r0 = r0(curr_idx, :);
    curr_r2 = r2(curr_idx, :);
    curr_len = sum(curr_idx);
    if curr_len > 100
        fprintf('!!!');
    end
    
    % find rotation between (r0 + r2) and (ray_in + ray_out)
    tmp_vec_a = curr_r0 + curr_r2;
    tmp_vec_a = bsxfun(@times, tmp_vec_a, 1./sqrt(sum(tmp_vec_a.^2, 2)));
    tmp_vec_b = ray_in + curr_ray_out;
    tmp_vec_b = tmp_vec_b / norm(tmp_vec_b);
    q1_axis = zeros(size(tmp_vec_a));
    for j = 1:curr_len
        q1_axis(j, :) = cross(tmp_vec_a(j, :), tmp_vec_b);
    end
    tmp_cos = tmp_vec_a * tmp_vec_b';
    q1_theta = atan2d(sqrt(sum(q1_axis.^2, 2)), tmp_cos);
    q1_axis = bsxfun(@times, q1_axis, 1./sqrt(sum(q1_axis.^2, 2)));
    q1 = [cosd(-q1_theta/2), bsxfun(@times, sind(-q1_theta/2), q1_axis)];
    
    % find rotation between rotated (r0 - r2) and (ray_in - ray_out)
    tmp_vec_a = quatrotate(q1, curr_r0 - curr_r2);
    tmp_vec_a = bsxfun(@times, tmp_vec_a, 1./sqrt(sum(tmp_vec_a.^2, 2)));
    tmp_vec_b = ray_in - curr_ray_out;
    tmp_vec_b = tmp_vec_b / norm(tmp_vec_b);
    q2_axis = zeros(size(tmp_vec_a));
    for j = 1:curr_len
        q2_axis(j, :) = cross(tmp_vec_a(j, :), tmp_vec_b);
    end
    tmp_cos = tmp_vec_a * tmp_vec_b';
    q2_theta = atan2d(sqrt(sum(q2_axis.^2, 2)), tmp_cos);
    q2_axis0 = (ray_in + curr_ray_out) / norm(ray_in + curr_ray_out);
    q2_axis = bsxfun(@times, q2_axis0, sign(q2_axis * q2_axis0'));
    q2 = [cosd(-q2_theta/2), bsxfun(@times, sind(-q2_theta/2), q2_axis)];
    
    % total quaternion
    q_all = quatmultiply(q1, q2);
    
    % crystal main axis
    main_axis = quatrotate(q_all, [0, 0, 1]);
    lambda = atan2d(main_axis(:, 2), main_axis(:, 1));
    zen = 90 - asind(main_axis(:, 3));
    w = exp(-(zen - crystal_zenith_mean).^2 / 2 / crystal_zenith_std^2) / crystal_zenith_std;
    halo_img(i) = sum(w);
end

%%
figure(2); clf;
scatter(r0_ll(:, 1), r0_ll(:, 2), 10, halo_img, 'filled');
axis equal; axis tight;
box on;

