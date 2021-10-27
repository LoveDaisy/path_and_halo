clear; close all; clc;

n = 1.31;
wedge_angle = 60;
norm_angle = 180 - wedge_angle;
wedge_n = [roty(-norm_angle / 2) * [0; 0; 1], roty(norm_angle / 2) * [0; 0; 1]]';

n_side = 2^6;
n_pix = nSide2nPix(n_side);

out_ray = nan(n_pix, 2);
quat_store = nan(n_pix, 4);

[ray_in_x, ray_in_y, ray_in_z] = pix2vec(n_side, 1:n_pix);
ray_in = -[ray_in_x', ray_in_y', ray_in_z'];
clear ray_in_x ray_in_y ray_in_z;

valid_entry_flag = ray_in * wedge_n(1, :)' < 0;
r1_store = refract(ray_in, wedge_n(1, :), 1, n);
valid_r1_flag = sum(abs(r1_store), 2) > 1e-4;
r2_store = refract(r1_store, wedge_n(2, :), n, 1);
valid_r2_flag = sum(abs(r2_store), 2) > 1e-4;

for i = 1:n_pix;
    if ~valid_entry_flag(i)
        out_ray(i, 2) = -1;
        continue;
    end

    r1 = r1_store(i, :);
    if norm(r1) < 1e-4
        out_ray(i, 2) = -.5;
        continue;
    end
    r2 = r2_store(i, :);
    if norm(r2) < 1e-4
        out_ray(i, 2) = 0;
        continue;
    end
    
    r1_axis = cross(r1, [1, 0, 0]);
    r1_theta = atan2d(norm(r1_axis), dot(r1, [1, 0, 0]));
    r1_axis = r1_axis / norm(r1_axis);
    r1_quat = [cosd(-r1_theta / 2), r1_axis * sind(-r1_theta / 2)];
    
    internal_r2 = quatrotate(r1_quat, r2);
    
    r2_axis = [1, 0, 0];
    r2_theta = atan2d(internal_r2(2), internal_r2(3));
    r2_quat = [cosd(-r2_theta / 2), r2_axis * sind(-r2_theta / 2)];
    
    final_quat = quatmultiply(r1_quat, r2_quat);
    new_r2 = quatrotate(final_quat, r2);
    
    quat_store(i, :) = final_quat;
    out_ray(i, :) = [atan2d(new_r2(3), new_r2(1)), 1];
end

%%
figure(1); clf;
subplot(2,1,1);
scatter(atan2d(ray_in(:, 2), ray_in(:, 1)), asind(ray_in(:, 3) ./ sqrt(sum(ray_in.^2, 2))), ...
    1, out_ray(:, 2), 'filled');
axis equal; axis tight;

subplot(2,1,2);
scatter(atan2d(ray_in(:, 2), ray_in(:, 1)), asind(ray_in(:, 3) ./ sqrt(sum(ray_in.^2, 2))), ...
    1, out_ray(:, 1), 'filled');
axis equal; axis tight;