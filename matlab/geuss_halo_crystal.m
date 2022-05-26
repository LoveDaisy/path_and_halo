clear; close all; clc;

test_img = im2double(imread('img_05.jpg'));
test_img = test_img(:, :, 2);
test_img = (test_img > 0.1) .* test_img;

test_img = imgaussfilt(test_img, 1);
test_img = imresize(test_img, 1/4, 'bilinear');
test_img = imgaussfilt(test_img, 1);
test_img = imresize(test_img, 1/2, 'bilinear');
test_img = test_img / max(test_img(:));
[img_hei, img_wid] = size(test_img);

k = img_hei / sqrt(2);
sun_ll = [0, 20];

axis_pdf = opt.make_axis_pdf();
crystal = opt.make_prism_crystal(1);

rp_list = {[1;2;3;5;1]; [1;2;5;3;1]};
rp_num = length(rp_list);
trace_list = cell(rp_num, 1);
config_list = cell(rp_num, 1);
fdf_list = cell(rp_num, 1);
for i = 1:length(rp_list)
    trace.fid = rp_list{i};
    trace_list{i} = trace;
    config_list{i} = opt.init_config(crystal, trace, sun_ll, 'GridLevel', 3);
    fdf_list{i} = @(rot) opt.crystal_system(rot, config_list{i}.ray_in_ll, crystal, trace_list{i});
end

%%
target_ll_store = nan(img_wid * img_hei, 3);
target_ll_cnt = 0;
rot_llr_store = nan(100 * img_wid * img_hei * rp_num, 4);
rot_llr_cnt = 0;
progress = utl.Progress(1 / (img_wid * img_hei), 0.02);
for x = 1:img_wid
    for y = 1:img_hei
        update_progress = progress.tik_and_show(1, 'process (%d,%d)', x, y);
        if test_img(y, x) < 1e-3
            continue;
        end
        az = atan2d(y - img_hei / 2, x - img_wid / 2) - 90;
        r = norm([x - img_wid / 2, y - img_hei / 2]);
        if r > img_hei / 2 - 0.5
            continue;
        end
        target_ll = [az + 180, asind(r / k) * 2 - 90];
        target_ll(1) = mod(target_ll(1) + 180, 360) - 180;
        target_ll(2) = mod(target_ll(2) + 90, 180) - 90;
        target_ll_store(target_ll_cnt + 1, 1:2) = target_ll;
        target_ll_store(target_ll_cnt + 1, 3) = test_img(y, x);
        target_ll_cnt = target_ll_cnt + 1;
        
        for ri = 1:rp_num
            rot_seeds = opt.find_cand_rot(config_list{ri}, target_ll, 'quat');
            if isempty(rot_seeds)
                continue;
            end
            
            [contour_store, rot_seeds, status] = opt.find_all_contours(target_ll, fdf_list{ri}, config_list{ri});
            if isempty(contour_store)
                continue;
            end
            
            for j = 1:length(contour_store)
                rot_contour = contour_store{j};
                % Compute weight
                [~, interp_cmp, interp_llr] = opt.compute_contour_weight(rot_contour, axis_pdf, config_list{ri});
                interp_w = interp_cmp(:, 4);
                
                interp_rot = geo.llr2quat(interp_llr);
                [d, ~, itp] = geo.distance_to_polyline(interp_rot, rot_seeds);
                
                seeds_cnt = length(d);
                rot_seeds_llr = geo.quat2llr(rot_seeds);
                rot_llr_store(rot_llr_cnt + (1:seeds_cnt), 1:3) = rot_seeds_llr;
                rot_llr_store(rot_llr_cnt + (1:seeds_cnt), 4) = test_img(y, x) * ...
                    interp_w(itp(:, 1)) .* exp(-d.^2 / 0.08^2);
                rot_llr_cnt = rot_llr_cnt + seeds_cnt;
            end
        end
    end
    
    tmp_data = rot_llr_store(1:rot_llr_cnt, :);
    tmp_data = repmat(tmp_data, 6, 1);
    for j = 1:6
        tmp_data(rot_llr_cnt * (j - 1) + (1:rot_llr_cnt), 3) = ...
            mod(tmp_data(rot_llr_cnt * (j - 1) + (1:rot_llr_cnt), 3) + 60 * j, 360);
    end
    tmp_data(:, 4) = tmp_data(:, 4) ./ max(cosd(tmp_data(:, 2)), 1e-6);
    figure(2); clf;
    utl.bubble_hist3(tmp_data, 10, 'scale', 2);
    title(sprintf('x: %d', x));
    axis equal;
    drawnow;
end
