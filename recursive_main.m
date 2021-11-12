clear; close all; clc;

crystal.face_norm = [0, 0, 1;     % face 1
            0, 0, -1;             % face 2
            -sqrt(3)/2, 1/2, 0;   % face 3
            0, 1, 0;              % face 4
            sqrt(3)/2, 1/2, 0;    % face 5
            sqrt(3)/2, -1/2, 0;   % face 6
            0, -1, 0;             % face 7
            -sqrt(3)/2, -1/2, 0]; % face 8
crystal.h = 10;
crystal.vtx = zeros(12, 3);
for i = 1:6
    crystal.vtx(i, :) = [cosd((4-i) * 60), sind((4-i)*60), crystal.h/2];
end
for i = 7:12
    crystal.vtx(i, :) = [cosd((4-i) * 60), sind((4-i)*60), -crystal.h/2];
end
crystal.face = {[6; 5; 4; 3; 2; 1];
                [7; 8; 9; 10; 11; 12];
                [1; 2; 8; 7];
                [2; 3; 9; 8];
                [3; 4; 10; 9];
                [4; 5; 11; 10];
                [5; 6; 12; 11];
                [6; 1; 7; 12]};
crystal.face_area = [3 * sqrt(3) / 2; 3 * sqrt(3) / 2; crystal.h * ones(6, 1)];

n = 1.31;
trace.n = [n; 1];
trace.fid = [3; 5];

crystal_zenith = [90, 0.2];  % mean, std
tmp_x = linspace(-90, 90, 50000);
tmp_pdf = exp(-(90 - tmp_x - crystal_zenith(1)).^2 / 2 / crystal_zenith(2)^2) / crystal_zenith(2);
tmp_total = (sum(tmp_pdf) * (tmp_x(2) - tmp_x(1)));
clear tmp_pdf tmp_x
axis_pdf = @(llq) (exp(-(90 - llq(:, 2) - crystal_zenith(1)).^2 / 2 / crystal_zenith(2)^2) / ...
    crystal_zenith(2) / tmp_total) * (1 / 360) * (1 / 360);

halo_vis_fun_helper = @(x, a, b) log10(x * b ./ (x + b) + a);
inv_vis_fun_helper = @(x, a, b) 1 ./ (1 ./ (10.^x - a) - 1/b);
halo_vis_fun = @(x) halo_vis_fun_helper(x, 2e-5, 1e-2);
inv_vis_fun = @(x) inv_vis_fun_helper(x, 2e-5, 1e-2);

%%
sun_altitude = 20;
sun_longitude = 0;
sun_ll = [sun_longitude, sun_altitude];
ray_in_xyz = ll2xyz_with_gradient(sun_ll);

config = init_config_3d(crystal, trace, sun_ll, 3);
line_color = colormap('lines');


%%
halo_img_res = 0.02;
halo_img_x = -60:halo_img_res:0;
halo_img_y = -30:halo_img_res:65;

halo_img = nan(length(halo_img_y), length(halo_img_x));
computed_img = zeros(size(halo_img));
level_img = zeros(size(halo_img));

checked_pix = 0;
progress_cnt = 0;
progress_bin = 0.0001;
update_progress = false;

recursive_level = max(min(floor(log2(5 / halo_img_res)), 7), 2);
% recursive_level = 5;
curr_step = 2^recursive_level;

recursive_stack = cell(32, 1);
next_stack = cell(32, 1);

next_stack_idx = 1;
next_stack{next_stack_idx}.box = [1, 1;
    1, 1 + curr_step;
    1 + curr_step, 1;
    1 + curr_step, 1 + curr_step];
next_stack{next_stack_idx}.level = recursive_level;
while next_stack_idx > 0
    curr_box = next_stack{next_stack_idx}.box;
    curr_level = next_stack{next_stack_idx}.level;
    curr_step = 2^curr_level;
    next_stack_idx = next_stack_idx - 1;
    
    recursive_stack_idx = 1;
    recursive_stack{recursive_stack_idx}.box = curr_box;
    recursive_stack{recursive_stack_idx}.level = curr_level;
    
    % Right box
    if curr_box(3, 1) + curr_step <= length(halo_img_x)
        x0 = curr_box(3, 1);
        y0 = curr_box(3, 2);
        next_box = [x0, y0;
            x0, y0 + curr_step;
            x0 + curr_step, y0;
            x0 + curr_step, y0 + curr_step];
        box_exist = false;
        for i = 1:next_stack_idx
            if sum(abs(next_stack{i}.box(:) - next_box(:))) < 1e-3
                box_exist = true;
                break;
            end
        end
        next_center = mean(next_box([1, 4], :));
        if ~box_exist && computed_img(next_center(2), next_center(1)) <= 0
            next_stack_idx = next_stack_idx + 1;
            next_stack{next_stack_idx}.box = next_box;
            next_stack{next_stack_idx}.level = curr_level;
        end
    elseif curr_box(3, 1) < length(halo_img_x)
        x0 = curr_box(3, 1);
        y0 = curr_box(3, 2);
        tmp_level = curr_level;
        while tmp_level > 0 && x0 + 2^tmp_level > length(halo_img_x)
            tmp_level = tmp_level - 1;
        end
        if tmp_level >= 0
            tmp_step = 2^tmp_level;
            next_box = [x0, y0;
                x0, y0 + tmp_step;
                x0 + tmp_step, y0;
                x0 + tmp_step, y0 + tmp_step];
            box_exist = false;
            for i = 1:next_stack_idx
                if sum(abs(next_stack{i}.box(:) - next_box(:))) < 1e-3
                    box_exist = true;
                    break;
                end
            end
            next_center = mean(next_box([1, 4], :));
            if ~box_exist && computed_img(next_center(2), next_center(1)) <= 0
                next_stack_idx = next_stack_idx + 1;
                next_stack{next_stack_idx}.box = next_box;
                next_stack{next_stack_idx}.level = tmp_level;
            end
        end
    end
    
    % Up box
    if curr_box(2, 2) + curr_step <= length(halo_img_y)
        x0 = curr_box(2, 1);
        y0 = curr_box(2, 2);
        next_box = [x0, y0;
            x0, y0 + curr_step;
            x0 + curr_step, y0;
            x0 + curr_step, y0 + curr_step];
        box_exist = false;
        for i = 1:next_stack_idx
            if sum(abs(next_stack{i}.box(:) - next_box(:))) < 1e-3
                box_exist = true;
                break;
            end
        end
        next_center = mean(next_box([1, 4], :));
        if ~box_exist && computed_img(next_center(2), next_center(1)) <= 0
            next_stack_idx = next_stack_idx + 1;
            next_stack{next_stack_idx}.box = next_box;
            next_stack{next_stack_idx}.level = curr_level;
        end
    elseif curr_box(2, 2) < length(halo_img_y)
        x0 = curr_box(2, 1);
        y0 = curr_box(2, 2);
        tmp_level = curr_level;
        while tmp_level > 0 && y0 + 2^tmp_level > length(halo_img_y)
            tmp_level = tmp_level - 1;
        end
        if tmp_level >= 0
            tmp_step = 2^tmp_level;
            next_box = [x0, y0;
                x0, y0 + tmp_step;
                x0 + tmp_step, y0;
                x0 + tmp_step, y0 + tmp_step];
            box_exist = false;
            for i = 1:next_stack_idx
                if sum(abs(next_stack{i}.box(:) - next_box(:))) < 1e-3
                    box_exist = true;
                    break;
                end
            end
            next_center = mean(next_box([1, 4], :));
            if ~box_exist && computed_img(next_center(2), next_center(1)) <= 0
                next_stack_idx = next_stack_idx + 1;
                next_stack{next_stack_idx}.box = next_box;
                next_stack{next_stack_idx}.level = tmp_level;
            end
        end
    end
    
    while recursive_stack_idx > 0
        curr_box = recursive_stack{recursive_stack_idx}.box;
        curr_level = recursive_stack{recursive_stack_idx}.level;
        curr_step = 2^curr_level;
        recursive_stack_idx = recursive_stack_idx - 1;
        curr_vtx_val = zeros(2, 2);
    
        for i = 1:4
            if curr_level < 1
                level_img(curr_box(i, 2), curr_box(i, 1)) = recursive_level - curr_level + 1;
            end
            rec_fun_call = 0;
            if isnan(halo_img(curr_box(i, 2), curr_box(i, 1)))
                lon = halo_img_x(curr_box(i, 1));
                lat = halo_img_y(curr_box(i, 2));
                [x_contour, y_val, jacobian] = find_axis_rot_contour(sun_ll, ...
                    [lon, lat], crystal, trace, 'config', config);
                rec_fun_call = rec_fun_call + 1;

                weight = 0;
                for k = 1:length(x_contour)
                    tmp_rot = x_contour{k};
                    tmp_j = jacobian{k};
                    [tmp_w, tmp_s, tmp_p, tmp_rot] = ...
                        compute_axis_rot_weight(tmp_rot, tmp_j, axis_pdf, crystal, trace, sun_ll);
                    weight = weight + tmp_w;
                end
                halo_img(curr_box(i, 2), curr_box(i, 1)) = weight;
                computed_img(curr_box(i, 2), curr_box(i, 1)) = 2;
                checked_pix = checked_pix + 1;
                progress_cnt = progress_cnt + 1 / numel(halo_img);
            end
            curr_vtx_val(i) = halo_img(curr_box(i, 2), curr_box(i, 1));
        end
        if rec_fun_call == 0 && curr_level < 1
            continue;
        end

        curr_vis_val = halo_vis_fun(curr_vtx_val);

        if curr_level > 1
            center_x = curr_box(1, 1) + curr_step / 2;
            center_y = curr_box(1, 2) + curr_step / 2;
            if isnan(halo_img(center_y, center_x))
                lon = halo_img_x(center_x);
                lat = halo_img_y(center_y);
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
                halo_img(center_y, center_x) = weight;
                computed_img(center_y, center_x) = 2;
%                 level_img(curr_box(i, 2), curr_box(i, 1)) = recursive_level - curr_level + 1;
                checked_pix = checked_pix + 1;
                progress_cnt = progress_cnt + 1 / numel(halo_img);
            end

            center_val = halo_vis_fun(halo_img(center_y, center_x));
        elseif curr_level > 0
            center_x = curr_box(1, 1) + curr_step / 2;
            center_y = curr_box(1, 2) + curr_step / 2;
            computed_img(center_y, center_x) = 1;
            center_val = nan;
        else
            center_val = mean(curr_vis_val(:));
        end

        flat_unit = 2e-5 / halo_img_res;
        diag_flat = abs(mean(curr_vis_val([2,3])) - mean(curr_vis_val([1,4]))) < 5 * flat_unit;
        center_flat = isnan(center_val) || abs(center_val - mean(curr_vis_val(:))) < 2 * flat_unit;
        if curr_level > 1 && ((diag_flat && center_flat && min(curr_vtx_val(:)) >= 1e-7) || ...
                (max(curr_vis_val(:)) - min(curr_vis_val(:)) < 1e-2))
            row_range = min(curr_box(:, 2)):max(curr_box(:, 2));
            col_range = min(curr_box(:, 1)):max(curr_box(:, 1));
            
            interp_img = inv_vis_fun(interp2(curr_vis_val, curr_level));
            origin_img = halo_img(row_range, col_range);
            
            pix_to_fill = isnan(origin_img);
            origin_img(pix_to_fill) = 0;
            curr_checked_pix = sum(pix_to_fill(:));
            halo_img(row_range, col_range) = interp_img .* pix_to_fill + origin_img .* ~pix_to_fill;
            level_img(row_range, col_range) = recursive_level - curr_level + 1;
            checked_pix = checked_pix + curr_checked_pix;
            progress_cnt = progress_cnt + curr_checked_pix / numel(halo_img);
        elseif curr_level > 0
            curr_step = curr_step / 2;
            curr_level = curr_level - 1;

            x0 = curr_box(1, 1);
            y0 = curr_box(1, 2);

            recursive_stack_idx = recursive_stack_idx + 1;
            recursive_stack{recursive_stack_idx}.box = [x0, y0;
                x0, y0 + curr_step;
                x0 + curr_step, y0;
                x0 + curr_step, y0 + curr_step];
            recursive_stack{recursive_stack_idx}.level = curr_level;

            recursive_stack_idx = recursive_stack_idx + 1;
            recursive_stack{recursive_stack_idx}.box = [x0, y0 + curr_step;
                x0, y0 + curr_step*2;
                x0 + curr_step, y0 + curr_step;
                x0 + curr_step, y0 + curr_step*2];
            recursive_stack{recursive_stack_idx}.level = curr_level;

            recursive_stack_idx = recursive_stack_idx + 1;
            recursive_stack{recursive_stack_idx}.box = [x0 + curr_step, y0;
                x0 + curr_step, y0 + curr_step;
                x0 + curr_step*2, y0;
                x0 + curr_step*2, y0 + curr_step];
            recursive_stack{recursive_stack_idx}.level = curr_level;

            recursive_stack_idx = recursive_stack_idx + 1;
            recursive_stack{recursive_stack_idx}.box = [x0 + curr_step, y0 + curr_step;
                x0 + curr_step, y0 + curr_step*2;
                x0 + curr_step*2, y0 + curr_step;
                x0 + curr_step*2, y0 + curr_step*2];
            recursive_stack{recursive_stack_idx}.level = curr_level;
        end
        
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
%             imagesc(halo_img_x, halo_img_y, computed_img);
            imagesc(halo_img_x, halo_img_y, level_img);
            axis equal; axis tight; axis xy;
            drawnow;
            
            progress_cnt = progress_cnt - floor(progress_cnt / progress_bin) * progress_bin;
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
