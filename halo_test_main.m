clear; close all; clc;

n = 1.31;
prism_h = 1;
face_norm = [0, 0, 1;     % face 1
    0, 0, -1;             % face 2
    -sqrt(3)/2, 1/2, 0;   % face 3
    0, 1, 0;              % face 4
    sqrt(3)/2, 1/2, 0;    % face 5
    sqrt(3)/2, -1/2, 0;   % face 6
    0, -1, 0;             % face 7
    -sqrt(3)/2, -1/2, 0]; % face 8
face_area = [sqrt(3) / 2; sqrt(3) / 2; prism_h; prism_h; prism_h; prism_h; prism_h; prism_h];

n_side = 2^6;
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

valid_idx = r0 * entry_norm' < 0;
r1(valid_idx, :) = refract(r0(valid_idx, :), entry_norm, 1, n);
valid_idx = valid_idx & (r1 * exit_norm' > 0);
r2(valid_idx, :) = refract(r1(valid_idx, :), exit_norm, n, 1);
valid_idx = valid_idx & sum(r2.^2, 2) > 1e-4;
r2(~valid_idx, :) = nan;

bending_angle = acosd(sum(r0 .* r2, 2));
bending_angle_max = max(bending_angle);
bending_angle_min = min(bending_angle);

%%
sun_altitude = 20;  % degree
sun_longitude = 180;
ray_in = -[cosd(sun_altitude) * cosd(sun_longitude), ...
    cosd(sun_altitude) * sind(sun_longitude), sind(sun_altitude)];
ray_out = -r0;

crystal_zenith_mean = 90;
crystal_zenith_std = 0.2;
dist_helper = @(zen) exp(-(zen - crystal_zenith_mean).^2 / ...
    2 / crystal_zenith_std^2) ./ crystal_zenith_std;
dist_x = linspace(-90, 90, 500);
dist_y = dist_helper(dist_x);
dist_total = sum(dist_y) * (dist_x(2) - dist_x(1));
crystal_zenith_dist = @(zen) dist_helper(zen) / dist_total;
clear dist_x dist_y;


img_dr = dr * 1;
img_wid = ceil(180 / img_dr) * 2;
img_hei = img_wid / 2;
lon_data = linspace(-180, 180, img_wid);
lat_data = linspace(-90, 90, img_hei);
[lon, lat] = meshgrid(lon_data, lat_data);
total_pix = numel(lon);
halo_img = zeros(img_hei, img_wid);
cm = repmat((0:255)', 1, 3) / 255;

progress_count = 0;
progress_bin = 0.002;
curr_len = 0;
search_size = zeros(total_pix, 1);
for i = 1:total_pix
    progress_count = progress_count + 1 / total_pix;
    if progress_count > progress_bin && curr_len > 0
        fprintf('processing %05.2f%%...\n', i / total_pix * 100);
        progress_count = progress_count - floor(progress_count / progress_bin) * progress_bin;
        figure(2); clf;
        imagesc(lon_data, lat_data, halo_img);
        colormap(cm);
        axis equal; axis tight;
        drawnow;
    end
    curr_ray_out = [cosd(lat(i)) * cosd(lon(i)), cosd(lat(i)) * sind(lon(i)), sind(lat(i))];
    curr_bending = acosd(dot(curr_ray_out, ray_in));
    if curr_bending < bending_angle_min || curr_bending > bending_angle_max
        curr_len = 0;
        continue;
    end
    
    curr_idx = abs(bending_angle - curr_bending) < img_dr & valid_idx;
    if ~any(curr_idx)
        curr_len = 0;
        continue;
    end
    
    curr_r0 = r0(curr_idx, :);
    curr_r2 = r2(curr_idx, :);
    curr_len = sum(curr_idx);
    search_size(i) = curr_len;
    
    % find rotation between r0 and ray_in
    tmp_vec_a = curr_r0;
    tmp_vec_b = ray_in;
    q1_axis = zeros(size(tmp_vec_a));
    for j = 1:curr_len
        q1_axis(j, :) = cross(tmp_vec_a(j, :), tmp_vec_b);
    end
    tmp_cos = tmp_vec_a * tmp_vec_b';
    q1_theta = atan2d(sqrt(sum(q1_axis.^2, 2)), tmp_cos);
    q1_axis = bsxfun(@times, q1_axis, 1./sqrt(sum(q1_axis.^2, 2)));
    q1 = [cosd(-q1_theta/2), bsxfun(@times, sind(-q1_theta/2), q1_axis)];
    
    % find rotation between rotated r2 and ray_out
    tmp_vec_a = quatrotate(q1, curr_r2);
    tmp_vec_a = tmp_vec_a - bsxfun(@times, tmp_vec_a * ray_in', ray_in);
    tmp_vec_a = bsxfun(@times, tmp_vec_a, 1./sqrt(sum(tmp_vec_a.^2, 2)));
    tmp_vec_b = curr_ray_out;
    tmp_vec_b = tmp_vec_b - (tmp_vec_b * ray_in') * ray_in;
    tmp_vec_b = tmp_vec_b / norm(tmp_vec_b);
    q2_axis = ray_in;
    q2_theta = acosd(tmp_vec_a * tmp_vec_b');
    q2 = [cosd(-q2_theta/2), bsxfun(@times, sind(-q2_theta/2), q2_axis)];
    q3 = [cosd(q2_theta/2), bsxfun(@times, sind(q2_theta/2), q2_axis)];
    curr_idx = quatrotate(q2, tmp_vec_a) * tmp_vec_b' > ...
        quatrotate(q3, tmp_vec_a) * tmp_vec_b';
    q2(~curr_idx, :) = q3(~curr_idx, :);
    
    % total quaternion
    q_all = quatmultiply(q1, q2);
    
    % projected face area
    proj_face_area = zeros(curr_len, size(face_norm, 1));
    for j = 1:size(face_norm, 1)
        proj_face_area(:, j) = -quatrotate(q_all, face_norm(j, :)) * ray_in';
    end
    geo_factor = max(proj_face_area(:, entry_face_idx), 0) ./ sum(max(proj_face_area, 0), 2);
    
    % crystal main axis
    main_axis = quatrotate(q_all, [0, 0, 1]);
    lambda = atan2d(main_axis(:, 2), main_axis(:, 1));
    zen = 90 - asind(main_axis(:, 3));
    w = crystal_zenith_dist(zen);
    halo_img(i) = sum(w .* geo_factor);
end

%%
figure(2); clf;
show_img = halo_img / max(halo_img(:));
imshow(show_img.^.85 * 2.5);
colormap(cm);
axis equal; axis tight;

