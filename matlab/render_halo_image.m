clear; close all; clc;

crystal = opt.make_prism_crystal(1);
% trace.fid = [1; 3; 2; 4; 5; 1];
trace.fid = [3; 5];

axis_pdf = opt.make_axis_pdf([0, 0.5]);
% axis_pdf = opt.make_axis_pdf();

sun_ll = [180, 15];
config = opt.init_config(crystal, trace, sun_ll);

fdf = @(rot) opt.crystal_system(rot, config.ray_in_ll, crystal, trace);

%%
tic;
halo_img = utl.generate_halo_image([22.5, 28.5], [-1, 1]-15, .1);

fun_eval_cnt = 0;
progress = utl.Progress(1 / (halo_img.x_length * halo_img.y_length), 0.005);

[xx, yy] = meshgrid(1:halo_img.x_length, 1:halo_img.y_length);
all_xy = [xx(:), yy(:)];
clear xx yy;
% all_xy = [21, 11];

for i = 1:size(all_xy, 1)
    w = all_xy(i, 1);
    h = all_xy(i, 2);
    target_ll = [halo_img.img_x(w), halo_img.img_y(h)];
    
    update_progress = progress.tik_and_show(1, 'process (%d,%d) = (%.2f,%.2f)', ...
        w, h, target_ll(1), target_ll(2));
    
    [contour_store, rot_seeds, status] = opt.find_all_contours(target_ll, fdf, config);
    if isempty(contour_store)
        continue;
    end
    
    contour_info = utl.ContourInfo;
    contour_info.add_candidate_seeds(rot_seeds);
    
    weight = 0;
    for j = 1:length(contour_store)
        rot_contour = contour_store{j};
        % Compute weight
        [curr_w, curr_cmp, curr_rot] = opt.compute_contour_weight(rot_contour, axis_pdf, config);
        weight = weight + curr_w;
        contour_info.add_contour(rot_contour, curr_cmp, curr_rot);
    end
    halo_img.img(h, w) = weight;
    
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
toc;

%%
figure(1); clf;
utl.show_halo_img(halo_img);
axis ij;
