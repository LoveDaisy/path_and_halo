clear; close all; clc;

crystal.face_norm = [0, 0, 1;     % face 1
            0, 0, -1;             % face 2
            -sqrt(3)/2, 1/2, 0;   % face 3
            0, 1, 0;              % face 4
            sqrt(3)/2, 1/2, 0;    % face 5
            sqrt(3)/2, -1/2, 0;   % face 6
            0, -1, 0;             % face 7
            -sqrt(3)/2, -1/2, 0]; % face 8
crystal.h = 1;
crystal.vtx = zeros(12, 3);
for i = 1:6
    crystal.vtx(i, :) = [cosd((4-i) * 60), sind((4-i)*60), crystal.h];
end
for i = 7:12
    crystal.vtx(i, :) = [cosd((4-i) * 60), sind((4-i)*60), -crystal.h];
end
crystal.face = {[6; 5; 4; 3; 2; 1];
                [7; 8; 9; 10; 11; 12];
                [1; 2; 8; 7];
                [2; 3; 9; 8];
                [3; 4; 10; 9];
                [4; 5; 11; 10];
                [5; 6; 12; 11];
                [6; 1; 7; 12]};
crystal.face_area = [3 * sqrt(3) / 2; 3 * sqrt(3) / 2; crystal.h * ones(6, 1) * 2];
crystal.n = 1.31;

trace.fid = [1; 3; 2; 4; 5; 1];

crystal_zenith = [90, 0.2];  % mean, std
tmp_x = linspace(-90, 90, 50000);
tmp_pdf = exp(-(90 - tmp_x - crystal_zenith(1)).^2 / 2 / crystal_zenith(2)^2) / crystal_zenith(2);
zen_total = (sum(tmp_pdf) * (tmp_x(2) - tmp_x(1)));
clear tmp_pdf tmp_x
axis_pdf = @(llr) (exp(-(90 - asind(sind(llr(:, 2))) - crystal_zenith(1)).^2 / 2 / crystal_zenith(2)^2) / ...
    crystal_zenith(2) / zen_total) * (1 / 360) * (1 / 360);

%%
halo_vis_fun_helper = @(x, a, b) log10(x * b ./ (x + b) + a);
inv_vis_fun_helper = @(x, a, b) 1 ./ (1 ./ (10.^x - a) - 1/b);
halo_vis_fun = @(x) halo_vis_fun_helper(x, 1e-5, 1e-2);
inv_vis_fun = @(x) inv_vis_fun_helper(x, 1e-5, 1e-2);

%%
sun_altitude = 10;
sun_longitude = 180;
sun_ll = [sun_longitude, sun_altitude];
ray_in_xyz = ll2xyz_with_gradient(sun_ll);

config = init_config_3d(crystal, trace, sun_ll, 3);
line_color = colormap('lines');
close all;

%%
halo_img_res = .5;
halo_img_x = 0:halo_img_res:10;
halo_img_y = -20:halo_img_res:20;

tic;
halo_img = nan(length(halo_img_y), length(halo_img_x));
checked_pix = 0;
progress_cnt = 0;
progress_bin = 0.001;
update_progress = false;
% for w = 1:length(halo_img_x)
%     for h = 1:length(halo_img_y)
for w = 1
    for h = 21
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
            [lon, lat], crystal, trace, 'config', config);
        if isempty(x_contour)
            continue;
        end
        
        weight = 0;
        interp_p_store = cell(size(x_contour));
        interp_rot_store = cell(size(x_contour));
        for k = 1:length(x_contour)
            curr_rot = x_contour{k};
            [tmp_w, tmp_s, tmp_p, tmp_rot] = ...
                compute_axis_rot_weight(curr_rot, axis_pdf, crystal, trace, sun_ll);
            interp_p_store{k} = [tmp_s, tmp_p];
            interp_rot_store{k} = tmp_rot;
            weight = weight + tmp_w;
        end
        halo_img(h, w) = weight;
        
%         if update_progress && weight > 1e-8
            figure(1); clf;
            f1_pos = get(gcf, 'position');
            imagesc(halo_img_x, halo_img_y, halo_vis_fun(halo_img));
            axis equal; axis tight; axis xy;
            title(sprintf('(%d,%d)\nw: %.4e', w, h, weight));
            drawnow;
            
            figure(2); clf;
            set(gcf, 'position', f1_pos + [f1_pos(3), 0, f1_pos(3), 0]);
            subplot(2,2,1);
            hold on;
            for k = 1:length(x_contour)
                plot(x_contour{k}(:, 2), x_contour{k}(:, 3), '-o');
                plot(interp_rot_store{k}(:, 2), interp_rot_store{k}(:, 3), '.-');
                title(sprintf('(%d,%d)\nw: %.4e', w, h, weight));
            end
            box on;
            lat_lim = get(gca, 'xlim'); roll_lim1 = get(gca, 'ylim');
            subplot(2,2,2);
            hold on;
            for k = 1:length(x_contour)
                plot(x_contour{k}(:, 1), x_contour{k}(:, 3), '-o');
                plot(interp_rot_store{k}(:, 1), interp_rot_store{k}(:, 3), '.-');
                title(sprintf('(%d,%d)\nw: %.4e', w, h, weight));
            end
            box on;
            lon_lim = get(gca, 'xlim'); roll_lim2 = get(gca, 'ylim');
            subplot(2,2,1);
            set(gca, 'ylim', [min(roll_lim1(1), roll_lim2(1)), max(roll_lim1(2), roll_lim2(2))]);
            subplot(2,2,2);
            set(gca, 'ylim', [min(roll_lim1(1), roll_lim2(1)), max(roll_lim1(2), roll_lim2(2))]);
            subplot(2,2,[3,4]);
            hold on;
            for k = 1:length(interp_p_store)
                plot(interp_p_store{k}(:, 1), interp_p_store{k}(:, 2), '.-k', 'linewidth', 2);
                plot(interp_p_store{k}(:, 1), interp_p_store{k}(:, 3), '.-', ...
                    'color', line_color(1, :));
                plot(interp_p_store{k}(:, 1), interp_p_store{k}(:, 4), '--', 'linewidth', 2, ...
                    'color', line_color(2, :));
                plot(interp_p_store{k}(:, 1), interp_p_store{k}(:, 5), ':', 'linewidth', 3, ...
                    'color', line_color(3, :));
                plot(interp_p_store{k}(:, 1), interp_p_store{k}(:, 6), '-', 'linewidth', 1.2, ...
                    'color', line_color(4, :));
            end
            plot([interp_p_store{1}(1,1), interp_p_store{1}(end,1)], [1, 1], ':k');
            set(gca, 'yscale', 'log', 'ylim', [1e-8, 1e4]);
            box on;
            drawnow;
%         end
    end
end
toc;

%%
figure(1); clf;
imagesc(halo_img_x, halo_img_y, halo_vis_fun(halo_img));
axis equal; axis tight; axis xy;

%%
% figure(3);
% imagesc([halo_img_x, -wrev(halo_img_x(1:end-1))], halo_img_y, ...
%     halo_vis_fun([halo_img, fliplr(halo_img(:, 1:end-1))]));
% axis xy; axis equal; axis tight;
