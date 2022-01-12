clear; close all; clc;

crystal = opt.make_prism_crystal(1);
% trace.fid = [1; 3; 2; 4; 5; 1];
trace.fid = [3; 5];

axis_pdf = geo.make_axis_pdf([1, 90, 0.2], []);

%%
sun_ll = [180, 10];
ray_in_ll = [sun_ll(1) + 180, -sun_ll(2)];
ray_in_xyz = geo.ll2xyz(ray_in_ll);

fdf = @(rot) opt.crystal_system(rot, ray_in_ll, crystal, trace);

config = opt.init_config(crystal, trace, sun_ll, 3);
line_color = colormap('lines');
close all;

%%
use_rot_quat = true;

tic;
halo_img = generate_halo_image([0, 10], [0, 30], 0.5);

fun_eval_cnt = 0;
checked_pix = 0;
progress_cnt = 0;
progress_bin = 0.005;
update_progress = false;
for w = 1:halo_img.x_length
    for h = 1:halo_img.y_length
% for w = 8
%     for h = 36
        ray_out_ll = [halo_img.img_x(w), halo_img.img_y(h)];

        checked_pix = checked_pix + 1;
        progress_cnt = progress_cnt + 1 / numel(halo_img.img);

        if progress_cnt > progress_bin
            fprintf('process (%d,%d) = (%.3f,%.3f), %05.2f%%\n', w, h, ray_out_ll(1), ray_out_ll(2), ...
                checked_pix / numel(halo_img.img) * 100);
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
        cmp_store = {};
        rot_store = {};
        contour_store = {};
        k = 1;
        while ~isempty(seed_rot)
            % Find contour
            [rot_contour, contour_status] = ode.find_contour(fdf, seed_rot(1, :), 'h', contour_h);
            seed_rot = seed_rot(2:end, :);
            fun_eval_cnt = fun_eval_cnt + contour_status.fun_eval_cnt;
            if isempty(rot_contour)
                continue;
            end

            if use_rot_quat
                contour_store{k} = geo.quat2llr(rot_contour);
            else
                contour_store{k} = rot_contour;
            end

            % Reduce seeds
            [seed_rot, reduce_status] = geo.reduce_pts_polyline(rot_contour, seed_rot, 'eps', reduce_eps);
            fun_eval_cnt = fun_eval_cnt + reduce_status.fun_eval_cnt;
            if size(rot_contour, 1) < 2
                continue;
            end

            % Compute weight
            [curr_w, curr_cmp, curr_rot] = opt.compute_contour_weight(rot_contour, axis_pdf, config);
            weight = weight + curr_w;
            cmp_store{k} = curr_cmp;
            rot_store{k} = curr_rot;
            k = k + 1;
        end
        halo_img.img(h, w) = weight;
        if k <= 1
            continue;
        end

        if update_progress && weight > 1e-8
            figure(1); clf;
            vis_range = [5e-4, 1];
            f1_pos = get(gcf, 'position');
            show_halo_img(halo_img, vis_range);
            title(sprintf('(%d,%d)\nw: %.4e', w, h, weight));
            drawnow;

            figure(2); clf;
            set(gcf, 'position', f1_pos + [f1_pos(3), 0, f1_pos(3), 0]);
            show_contour_info(contour_store, rot_store, cmp_store);
            drawnow;
        end
    end
end
toc;

%%
figure(1); clf;
vis_range = [5e-4, 1];
show_halo_img(halo_img, vis_range);

