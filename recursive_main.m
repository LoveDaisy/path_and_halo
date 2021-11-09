% clear; close all; clc;
% 
% crystal.face_norm = [0, 0, 1;     % face 1
%             0, 0, -1;             % face 2
%             -sqrt(3)/2, 1/2, 0;   % face 3
%             0, 1, 0;              % face 4
%             sqrt(3)/2, 1/2, 0;    % face 5
%             sqrt(3)/2, -1/2, 0;   % face 6
%             0, -1, 0;             % face 7
%             -sqrt(3)/2, -1/2, 0]; % face 8
% crystal.h = 10;
% crystal.vtx = zeros(12, 3);
% for i = 1:6
%     crystal.vtx(i, :) = [cosd((4-i) * 60), sind((4-i)*60), crystal.h/2];
% end
% for i = 7:12
%     crystal.vtx(i, :) = [cosd((4-i) * 60), sind((4-i)*60), -crystal.h/2];
% end
% crystal.face = {[6; 5; 4; 3; 2; 1];
%                 [7; 8; 9; 10; 11; 12];
%                 [1; 2; 8; 7];
%                 [2; 3; 9; 8];
%                 [3; 4; 10; 9];
%                 [4; 5; 11; 10];
%                 [5; 6; 12; 11];
%                 [6; 1; 7; 12]};
% crystal.face_area = [3 * sqrt(3) / 2; 3 * sqrt(3) / 2; crystal.h * ones(6, 1)];
% 
% n = 1.31;
% trace.n = [n; 1];
% trace.fid = [3; 5];
% 
% crystal_zenith = [90, 0.2];  % mean, std
% tmp_x = linspace(-90, 90, 50000);
% tmp_pdf = exp(-(90 - tmp_x - crystal_zenith(1)).^2 / 2 / crystal_zenith(2)^2) / crystal_zenith(2);
% tmp_total = (sum(tmp_pdf) * (tmp_x(2) - tmp_x(1)));
% clear tmp_pdf tmp_x
% axis_pdf = @(llq) (exp(-(90 - llq(:, 2) - crystal_zenith(1)).^2 / 2 / crystal_zenith(2)^2) / ...
%     crystal_zenith(2) / tmp_total) * (1 / 360) * (1 / 360);
% 
% halo_vis_fun_helper = @(x, a, b) log10(x * b ./ (x + b) + a);
% halo_vis_fun = @(x) halo_vis_fun_helper(x, 2e-5, 1e-2);
% 
% %%
% sun_altitude = 20;
% sun_longitude = 0;
% sun_ll = [sun_longitude, sun_altitude];
% ray_in_xyz = ll2xyz_with_gradient(sun_ll);
% 
% config = init_config_3d(crystal, trace, sun_ll, 3);
% line_color = colormap('lines');


%%
halo_img_x = linspace(-40,0,257);
halo_img_y = linspace(-25,55,513);
halo_img = nan(length(halo_img_y), length(halo_img_x));
computed_img = false(size(halo_img));
checked_pix = 0;
progress_cnt = 0;
progress_bin = 0.001;
update_progress = false;

recursive_level = 5;
max_step = 2^recursive_level;
for wi = 1:max_step:length(halo_img_x)-1
    for hi = 1:max_step:length(halo_img_y)-1
        fprintf('(%d,%d)\n', wi, hi);
        curr_step = max_step;
        
        recursive_stack = cell(32, 1);
        stack_idx = 1;
        recursive_stack{stack_idx}.box = [wi, hi;
            wi, hi + curr_step;
            wi + curr_step, hi;
            wi + curr_step, hi + curr_step];
        recursive_stack{stack_idx}.level = recursive_level;
        while stack_idx > 0
            curr_box = recursive_stack{stack_idx}.box;
            curr_level = recursive_stack{stack_idx}.level;
            curr_step = 2^curr_level;
            stack_idx = stack_idx - 1;
            curr_val = zeros(2, 2);
            
            if progress_cnt > progress_bin
                for i = 1:recursive_level - curr_level
                    fprintf('-');
                end
                fprintf('(%d,%d): %05.2f%%\n', curr_box(1,1), curr_box(1,2), ...
                    checked_pix / numel(halo_img) * 100);
                figure(3); clf;
                subplot(1,2,1)
                imagesc(halo_img_x, halo_img_y, halo_vis_fun(halo_img));
                axis equal; axis tight; axis xy;
                subplot(1,2,2)
                imagesc(halo_img_x, halo_img_y, computed_img);
                axis equal; axis tight; axis xy;
                drawnow;
                
                progress_cnt = progress_cnt - floor(progress_cnt / progress_bin) * progress_bin;
            end
            
            for i = 1:4
                if isnan(halo_img(curr_box(i, 2), curr_box(i, 1)))
                    lon = halo_img_x(curr_box(i, 1));
                    lat = halo_img_y(curr_box(i, 2));
                    [x_contour, y_val, jacobian] = find_axis_rot_contour(sun_ll, ...
                        [lon, lat], crystal, trace, 'config', config);
                    
                    weight = 0;
                    for k = 1:length(x_contour)
                        tmp_rot = x_contour{k};
                        tmp_j = jacobian{k};
                        [tmp_w, tmp_s, tmp_p, tmp_rot] = ...
                            compute_axis_rot_weight(tmp_rot, tmp_j, axis_pdf, crystal, trace, sun_ll);
                        weight = weight + tmp_w;
                    end
                    halo_img(curr_box(i, 2), curr_box(i, 1)) = weight;
                    computed_img(curr_box(i, 2), curr_box(i, 1)) = true;
                    checked_pix = checked_pix + 1;
                    progress_cnt = progress_cnt + 1 / numel(halo_img);
                end
                curr_val(i) = halo_img(curr_box(i, 2), curr_box(i, 1));
            end
            curr_vis_val = halo_vis_fun(curr_val);
            curr_n = cross([curr_step, 0, curr_vis_val(1, 2)] - [0, 0, curr_vis_val(1, 1)], ...
                [0, curr_step, curr_vis_val(2, 1)] - [0, 0, curr_vis_val(1, 1)]);
            curr_n = curr_n / norm(curr_n);
            vd = [curr_step, curr_step, curr_vis_val(2, 2)] - [0, 0, curr_vis_val(1, 1)];
            err_sin_q = dot(vd, curr_n) / norm(vd);
            if (abs(err_sin_q) < 1e-3 && min(curr_val(:)) >= 1e-6) || ...
                    (max(curr_vis_val(:)) - min(curr_vis_val(:)) < 1e-2)
                row_range = min(curr_box(:, 2)):max(curr_box(:, 2));
                col_range = min(curr_box(:, 1)):max(curr_box(:, 1));
                interp_img = interp2(curr_val, curr_level);
                origin_img = halo_img(row_range, col_range);
                pix_to_fill = isnan(origin_img);
                origin_img(pix_to_fill) = 0;
                curr_checked_pix = sum(pix_to_fill(:));
                halo_img(row_range, col_range) = interp_img .* pix_to_fill + origin_img .* ~pix_to_fill;
                checked_pix = checked_pix + curr_checked_pix;
                progress_cnt = progress_cnt + curr_checked_pix / numel(halo_img);
            elseif curr_level > 0
                curr_step = curr_step / 2;
                curr_level = curr_level - 1;
                
                box_x0 = curr_box(1, 1);
                box_y0 = curr_box(1, 2);
                
                stack_idx = stack_idx + 1;
                recursive_stack{stack_idx}.box = [box_x0, box_y0;
                    box_x0, box_y0 + curr_step;
                    box_x0 + curr_step, box_y0;
                    box_x0 + curr_step, box_y0 + curr_step];
                recursive_stack{stack_idx}.level = curr_level;
                
                stack_idx = stack_idx + 1;
                recursive_stack{stack_idx}.box = [box_x0, box_y0 + curr_step;
                    box_x0, box_y0 + curr_step*2;
                    box_x0 + curr_step, box_y0 + curr_step;
                    box_x0 + curr_step, box_y0 + curr_step*2];
                recursive_stack{stack_idx}.level = curr_level;
                
                stack_idx = stack_idx + 1;
                recursive_stack{stack_idx}.box = [box_x0 + curr_step, box_y0;
                    box_x0 + curr_step, box_y0 + curr_step;
                    box_x0 + curr_step*2, box_y0;
                    box_x0 + curr_step*2, box_y0 + curr_step];
                recursive_stack{stack_idx}.level = curr_level;
                
                stack_idx = stack_idx + 1;
                recursive_stack{stack_idx}.box = [box_x0 + curr_step, box_y0 + curr_step;
                    box_x0 + curr_step, box_y0 + curr_step*2;
                    box_x0 + curr_step*2, box_y0 + curr_step;
                    box_x0 + curr_step*2, box_y0 + curr_step*2];
                recursive_stack{stack_idx}.level = curr_level;
            else
                continue;
            end
        end
    end
end


%%
figure(1); clf;
imagesc(halo_img_x, halo_img_y, halo_vis_fun(halo_img));
axis equal; axis tight; axis xy;

%%
figure(2);
imagesc([halo_img_x, -wrev(halo_img_x(1:end-1))], halo_img_y, ...
    halo_vis_fun([halo_img, fliplr(halo_img(:, 1:end-1))]));
axis xy; axis equal; axis tight;
