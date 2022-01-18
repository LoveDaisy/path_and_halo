clear; close all; clc;

crystal = opt.make_prism_crystal(1);
trace.fid = [1; 3; 2; 4; 5; 1];
% trace.fid = [3; 5];

axis_pdf = geo.make_axis_pdf([1, 90, 0.2], []);
% axis_pdf = geo.make_axis_pdf([0, 0, 0], []);

sun_ll = [0, 15];
ray_in_ll = [sun_ll(1) + 180, -sun_ll(2)];
ray_in_xyz = geo.ll2xyz(ray_in_ll);

fdf = @(rot) opt.crystal_system(rot, ray_in_ll, crystal, trace);

config = opt.init_config(crystal, trace, sun_ll, 3);
line_color = colormap('lines');
close all;

%%
use_rot_quat = true;

tic;
halo_img = utl.generate_halo_image([-1.4, 1.4], [12, 18], .1);

fun_eval_cnt = 0;
progress = utl.Progress(1 / (halo_img.x_length * halo_img.y_length), 0.005);
for w = 1:halo_img.x_length
    for h = 1:halo_img.y_length
% for w = 17
%     for h = 29
        ray_out_ll = [halo_img.img_x(w), halo_img.img_y(h)];

        update_progress = progress.tik_and_show(1, 'process (%d,%d) = (%.2f,%.2f)', ...
            w, h, ray_out_ll(1), ray_out_ll(2));

        % Find seed rotation
        if use_rot_quat
            cand_rot = opt.find_cand_rot(config, ray_out_ll, 'quat');
            contour_h = 0.05;
            reduce_eps = 0.05;
        else
            cand_rot = opt.find_cand_rot(config, ray_out_ll, 'llr');
            contour_h = 1;
            reduce_eps = 0.5;
        end
        if isempty(cand_rot)
            continue;
        end

        weight = 0;
        contour_info = utl.ContourInfo;
        contour_info.add_candidate_seeds(cand_rot);
        while ~isempty(cand_rot)
            % Find initial solution
            [init_rot, sol_status] = ode.find_solution(fdf, cand_rot(1, :), [ray_out_ll, 1], 'eps', 1e-8);
            fun_eval_cnt = fun_eval_cnt + sol_status.fun_eval_cnt;
            cand_rot = cand_rot(2:end, :);
            if ~sol_status.finish || (use_rot_quat && init_rot(1) < 0)
                continue;
            end
            
            % Check if it locates on previous contour
            duplicated = false;
            for i = 1:contour_info.total_cnt
                [init_rot, dup_status] = geo.reduce_pts_polyline(contour_info.contour_store{i}, init_rot, ...
                    'eps', reduce_eps);
                fun_eval_cnt = fun_eval_cnt + dup_status.fun_eval_cnt;
                if isempty(init_rot)
                    duplicated = true;
                    break;
                end
            end
            if duplicated
                continue;
            end

            % Find contour
            [rot_contour, contour_status] = ode.find_contour(fdf, init_rot, 'h', contour_h);
            fun_eval_cnt = fun_eval_cnt + contour_status.fun_eval_cnt;

            if isempty(rot_contour)
                continue;
            end

            % Reduce seeds
            [cand_rot, reduce_status] = geo.reduce_pts_polyline(rot_contour, cand_rot, ...
                'eps', config.dr * 2.5, 'jac_fun', fdf);
            fun_eval_cnt = fun_eval_cnt + reduce_status.fun_eval_cnt;

            % Compute weight
            [curr_w, curr_cmp, curr_rot] = opt.compute_contour_weight(rot_contour, axis_pdf, config);
            weight = weight + curr_w;
            contour_info.add_contour(rot_contour, curr_cmp, curr_rot);
        end
        halo_img.img(h, w) = weight;
        if contour_info.total_cnt <= 1
            continue;
        end

        if update_progress
            figure(1); clf;
            utl.show_halo_img(halo_img);
            title(sprintf('(%d,%d)\nw: %.4e', w, h, weight));
            axis ij;
            drawnow;

            figure(2); clf;
            contour_info.display_info();
            drawnow;
        end
    end
end
toc;

%%
figure(1); clf;
utl.show_halo_img(halo_img);
axis ij;

