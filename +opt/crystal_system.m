function [ray_out_ll, g_rot] = crystal_system(rot_llr, ray_in_ll, crystal, trace)
% INPUT
%   rot_llr:            n*3
%   ray_in_ll:          1*2
%   crystal:            struct
%       .face_norm      m*3
%   trace:              struct
%       .fid            k*1
%       .n:             k*1

num = size(rot_llr, 1);

xyz0 = geo.ll2xyz(ray_in_ll);
q1 = rot_llr(:, 1);
q2 = rot_llr(:, 2);
q3 = rot_llr(:, 3);
cq1 = cosd(q1);
sq1 = sind(q1);
cq2 = cosd(q2);
sq2 = sind(q2);
cq3 = cosd(q3);
sq3 = sind(q3);

g_r_mat = zeros(3, 3, num, 3);
g_r_mat(:, :, :, 1) = reshape( ...
    [-cq1 .* cq3 + sq1 .* sq2 .* sq3, -sq1 .* cq3 - cq1 .* sq2 .* sq3, zeros(num, 1), ...
    sq1 .* sq2 .* cq3 + cq1 .* sq3, -cq1 .* sq2 .* cq3 + sq1 .* sq3, zeros(num, 1), ...
    -sq1 .* cq2, cq1 .* cq2, zeros(num, 1)]', [3, 3, num, 1]) * pi / 180;
g_r_mat(:, :, :, 2) = reshape( ...
    [-cq1 .* cq2 .* sq3, -sq1 .* cq2 .* sq3, -sq2 .* sq3, ...
    -cq1 .* cq2 .* cq3, -sq1 .* cq2 .* cq3, -sq2 .* cq3, ...
    -cq1 .* sq2, -sq1 .* sq2, cq2]', [3, 3, num, 1]) * pi / 180;
g_r_mat(:, :, :, 3) = reshape( ...
    [-cq1 .* sq2 .* cq3 + sq1 .* sq3, -sq1 .* sq2 .* cq3 - cq1 .* sq3, cq2 .* cq3, ...
    sq1 .* cq3 + cq1 .* sq2 .* sq3, -cq1 .* cq3 + sq1 .* sq2 .* sq3, -cq2 .* sq3, ...
    zeros(num, 3)]', [3, 3, num, 1]) * pi / 180;

ray_in_xyz = zeros(num, 3);
rot_mat_store = zeros(3, 3, num);
for i = 1:num
    curr_rot_mat = rotz(90 + rot_llr(i, 1)) * rotx(90 - rot_llr(i, 2)) * rotz(rot_llr(i, 3));
    rot_mat_store(:, :, i) = curr_rot_mat;
    ray_in_xyz(i, :) = xyz0 * curr_rot_mat;
end

[ray_out_xyz, g_out] = opt.trace_ray_xyz(ray_in_xyz, crystal, trace);

g_xyz = zeros(3, 3, num);
for i = 1:num
    g_xyz1 = [xyz0 * g_r_mat(:, :, i, 1); xyz0 * g_r_mat(:, :, i, 2); xyz0 * g_r_mat(:, :, i, 3)]';
    g_xyz(:, :, i) = [ray_out_xyz(i, :) * g_r_mat(:, :, i, 1)'; ray_out_xyz(i, :) * g_r_mat(:, :, i, 2)'; ...
        ray_out_xyz(i, :) * g_r_mat(:, :, i, 3)']' +  rot_mat_store(:, :, i) * g_out(:, :, i) * g_xyz1;
    ray_out_xyz(i, :) = ray_out_xyz(i, :) * rot_mat_store(:, :, i)';
end

[ray_out_ll, g_ll] = geo.xyz2ll(ray_out_xyz);
g_rot = zeros(2, 3, num);
for i = 1:num
    g_rot(:, :, i) = g_ll(:, :, i) * g_xyz(:, :, i);
end
end