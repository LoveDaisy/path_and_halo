clear; close all; clc;

face_norm = [0, 0, 1;     % face 1
    0, 0, -1;             % face 2
    -sqrt(3)/2, 1/2, 0;   % face 3
    0, 1, 0;              % face 4
    sqrt(3)/2, 1/2, 0;    % face 5
    sqrt(3)/2, -1/2, 0;   % face 6
    0, -1, 0;             % face 7
    -sqrt(3)/2, -1/2, 0]; % face 8
crystal_h = 10;
face_area = [3 * sqrt(3) / 2; 3 * sqrt(3) / 2; crystal_h * ones(6, 1)];
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
axis_pdf = @(llq) (exp(-(90 - llq(:, 2) - crystal_zenith(1)).^2 / 2 / crystal_zenith(2)^2) / ...
    crystal_zenith(2) / tmp_total) * (1 / 360) * (1 / 360);

%%
sun_altitude = 20;
sun_longitude = 0;
sun_ll = [sun_longitude, sun_altitude];
ray_in_xyz = ll2xyz_with_gradient(sun_ll);

config = init_config_3d(crystal_face_norm, crystal_n, sun_ll, 3);

%%
halo_img_x = -4:.01:0;
halo_img_y = 41.5:.01:42.5;
halo_img = zeros(length(halo_img_y), length(halo_img_x));
checked_pix = 0;
progress_cnt = 0;
progress_bin = 0.001;
update_progress = false;
for w = 53:length(halo_img_x)
    for h = 1:length(halo_img_y)
% for w = 315
%     for h = 36
        lon = halo_img_x(w);
        lat = halo_img_y(h);

        ray_out_xyz = ll2xyz_with_gradient([lon, lat]);
        curr_bending_angle = acosd(dot(ray_out_xyz, ray_in_xyz));
        curr_target_input_output = [ray_in_xyz, ray_out_xyz];
        
        checked_pix = checked_pix + 1;
        progress_cnt = progress_cnt + 1 / numel(halo_img);
        
        if progress_cnt > progress_bin
            fprintf('process (%d,%d) = (%.3f,%.3f), %05.2f%%\n', w, h, lon, lat, ...
                checked_pix / numel(halo_img) * 100);
            progress_cnt = progress_cnt - floor(progress_cnt / progress_bin) * progress_bin;
            update_progress = true;
        else
            update_progress = false;
        end

        [x_contour, y_val, jacobian] = find_axis_rot_contour(sun_ll, ...
            [lon, lat], crystal_face_norm, crystal_n, 'config', config);
        if isempty(x_contour)
            continue;
        end
        
        weight = 0;
        interp_p_store = cell(size(x_contour));
        interp_rot_store = cell(size(x_contour));
        for k = 1:length(x_contour)
            curr_rot = x_contour{k};
            curr_j = jacobian{k};
            [tmp_w, tmp_s, tmp_p, tmp_rot] = ...
                compute_axis_rot_weight(curr_rot, curr_j, axis_pdf, face_norm, face_area, ...
                entry_face_idx, sun_ll);
            interp_p_store{k} = [tmp_s, tmp_p];
            interp_rot_store{k} = tmp_rot;
            weight = weight + tmp_w;
        end
        halo_img(h, w) = weight;
        
        if weight > 1e-8 && update_progress
            figure(1); clf;
            f1_pos = get(gcf, 'position');
            imagesc(halo_img_x, halo_img_y, log10(halo_img * 1e-2 ./ (halo_img + 1e-2) + 5e-5));
            axis equal; axis tight; axis xy;
            title(sprintf('(%d,%d)\nw: %.4e', w, h, weight));
            drawnow;
            
            figure(2); clf;
            set(gcf, 'position', f1_pos + [f1_pos(3), 0, 0, 0]);
            subplot(2,1,1);
            hold on;
            for k = 1:length(x_contour)
                plot(x_contour{k}(:, 2), x_contour{k}(:, 3), '-o');
                plot(interp_rot_store{k}(:, 2), interp_rot_store{k}(:, 3), '.-');
                title(sprintf('(%d,%d)\nw: %.4e', w, h, weight));
            end
            subplot(2,1,2);
            hold on;
            for k = 1:length(interp_p_store)
                plot(interp_p_store{k}(:, 1), interp_p_store{k}(:, 2), '.-');
            end
            drawnow;
        end
    end
end

%%
figure(1); clf;
imagesc(halo_img_x, halo_img_y, log10(halo_img * 1e-2 ./ (halo_img + 1e-2) + 5e-5));
axis equal; axis tight; axis xy;
drawnow;

