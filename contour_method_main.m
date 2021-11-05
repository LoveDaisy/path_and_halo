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

entry_face_idx = 3;
exit_face_idx = 5;

crystal_face_norm = face_norm([entry_face_idx, exit_face_idx], :);
crystal_n = [n, 1];

crystal_zenith = [90, 0.2];  % mean, std
tmp_x = linspace(-90, 90, 50000);
tmp_pdf = exp(-(90 - tmp_x - crystal_zenith(1)).^2 / 2 / crystal_zenith(2)^2) / crystal_zenith(2);
tmp_total = (sum(tmp_pdf) * (tmp_x(2) - tmp_x(1)));
clear tmp_pdf tmp_x
axis_pdf = @(lon, lat) (exp(-(90 - lat - crystal_zenith(1)).^2 / 2 / crystal_zenith(2)^2) / ...
    crystal_zenith(2) / tmp_total) * (1 / 360) * (1 / 360);

%%
sun_altitude = 20;
sun_longitude = 0;
sun_ll = [sun_longitude, sun_altitude];
ray_in_xyz = ll2xyz_with_gradient(sun_ll);

config = init_config_3d(crystal_face_norm, crystal_n, sun_ll, 3);

halo_img_x = -10:.5:10;
halo_img_y = -30:.5:0;
halo_img = zeros(length(halo_img_y), length(halo_img_x));
checked_pix = 0;
progress_cnt = 0;
progress_bin = 0.001;
% for w = 1:length(halo_img_x)
%     for h = 1:length(halo_img_y)
% for w = 240:256
%     for h = 90:150
for w = 20
    for h = 40
        lon = halo_img_x(w);
        lat = halo_img_y(h);
        ray_out_xyz = ll2xyz_with_gradient([lon, lat]);
        curr_bending_angle = acosd(dot(ray_out_xyz, ray_in_xyz));
        curr_target_input_output = [ray_in_xyz, ray_out_xyz];
        
        checked_pix = checked_pix + 1;
        progress_cnt = progress_cnt + 1 / numel(halo_img);

        [x_contour, y_val, jacobian] = find_axis_rot_contour(sun_ll, ...
            [lon, lat], crystal_face_norm, crystal_n, 'config', config);
        if isempty(x_contour)
            continue;
        end
    end
end

