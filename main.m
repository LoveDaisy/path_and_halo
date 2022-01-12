clear; close all; clc;

crystal = opt.make_prism_crystal(1);
% trace.fid = [1; 3; 2; 4; 5; 1];
trace.fid = [3; 5];

crystal_zenith = [90, 0.2]; % mean, std
tmp_x = linspace(-90, 90, 50000);
tmp_pdf = exp(- (90 - tmp_x - crystal_zenith(1)).^2/2 / crystal_zenith(2)^2) / crystal_zenith(2);
zen_total = (sum(tmp_pdf) * (tmp_x(2) - tmp_x(1)));
clear tmp_pdf tmp_x
axis_pdf = @(llr) (exp(- (90 - asind(sind(llr(:, 2))) - crystal_zenith(1)).^2/2 / crystal_zenith(2)^2) / ...
    crystal_zenith(2) / zen_total) * (1/360) * (1/360);

halo_vis_fun_helper = @(x, a, b) log10(x * b ./ (x + b) + a);
inv_vis_fun_helper = @(x, a, b) 1 ./ (1 ./ (10.^x - a) - 1 / b);
halo_vis_fun = @(x) halo_vis_fun_helper(x, 1e-5, 1e-2);
inv_vis_fun = @(x) inv_vis_fun_helper(x, 1e-5, 1e-2);

%%
sun_altitude = 10;
sun_longitude = 180;
sun_ll = [sun_longitude, sun_altitude];
ray_in_ll = [sun_longitude + 180, -sun_altitude];
ray_in_xyz = geo.ll2xyz(ray_in_ll);

fdf = @(rot) opt.crystal_system(rot, ray_in_ll, crystal, trace);

config = opt.init_config(crystal, trace, sun_ll, 3);
line_color = colormap('lines');
close all;

%%
halo_img_res = .5;
halo_img_x = 0:halo_img_res:10;
halo_img_y = 0:halo_img_res:30;

use_rot_quat = true;

tic;
halo_img = nan(length(halo_img_y), length(halo_img_x));
fun_eval_cnt = 0;
checked_pix = 0;
progress_cnt = 0;
progress_bin = 0.005;
update_progress = false;
for w = 1:length(halo_img_x)
    for h = 1:length(halo_img_y)
% for w = 8
%     for h = 36
        lon = halo_img_x(w);
        lat = halo_img_y(h);
        ray_out_ll = [lon, lat];

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

        % Find seed rotation
        if use_rot_quat
            [seed_rot, seed_status] = opt.find_seed_rot(config, ray_out_ll, 'quat');
            contour_h = 0.05;
            reduce_eps = 0.05;
        else
            [seed_rot, seed_status] = opt.find_seed_rot(config, ray_out_ll, 'llr');
            contour_h = 1;
            reduce_eps = 0.5;
        end
        fun_eval_cnt = fun_eval_cnt + seed_status.fun_eval_cnt;
        if isempty(seed_rot) || size(seed_rot, 1) > 100
            continue;
        end
        
        weight = 0;
        interp_p_store = {};
        interp_rot_store = {};
        x_contour = {};
        k = 1;
        while ~isempty(seed_rot)
            % Find contour
            [rot_contour, contour_status] = ode.find_contour_fdf(fdf, seed_rot(1, :), 'h', contour_h);
            seed_rot = seed_rot(2:end, :);
            fun_eval_cnt = fun_eval_cnt + contour_status.fun_eval_cnt;
            if isempty(rot_contour)
                continue;
            end

            if use_rot_quat
                x_contour{k} = geo.quat2llr(rot_contour);
            else
                x_contour{k} = rot_contour;
            end
            
            % Reduce seeds
            [seed_rot, reduce_status] = geo.reduce_pts_polyline(rot_contour, seed_rot, 'eps', reduce_eps);
            fun_eval_cnt = fun_eval_cnt + reduce_status.fun_eval_cnt;
            if size(rot_contour, 1) < 2
                continue;
            end
            
            [curr_w, curr_cmp, curr_rot] = opt.compute_contour_weight(rot_contour, axis_pdf, config);
            weight = weight + curr_w;
            interp_p_store{k} = curr_cmp;
            interp_rot_store{k} = curr_rot;
            k = k + 1;
        end
        halo_img(h, w) = weight;
        if isempty(x_contour)
            continue;
        end
        
        if update_progress && weight > 1e-8
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
        end
    end
end
toc;

%%
% figure(1); clf;
% imagesc(halo_img_x, halo_img_y, halo_vis_fun(halo_img));
% axis equal; axis tight; axis xy;

%%
% figure(3);
% imagesc([halo_img_x, -wrev(halo_img_x(1:end-1))], halo_img_y, ...
%     halo_vis_fun([halo_img, fliplr(halo_img(:, 1:end-1))]));
% axis xy; axis equal; axis tight;
