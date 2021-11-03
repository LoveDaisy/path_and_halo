function [ray_out_ll, g_ll_rot] = crystal_system_with_gradient(rot_llq, ray_in_ll, face_norm, n)
% INPUT
%   rot_llq:            n*3
%   ray_in_ll:          1*2
%   face_norm:          k*3
%   n:                  k*1

num = size(ray_in_ll, 1);

xyz0 = ll2xyz_with_gradient(ray_in_ll);

ray_in_xyz = zeros(num, 3);
for i = 1:num
    curr_rot_mat = rotz(90 + rot_llq(1)) * rotx(90 - rot_llq(2)) * rotz(rot_llq(3));
    ray_in_xyz(i, :) = xyz0 * curr_rot_mat;
end

[ray_out_xyz, ~, g_out, ~] = trace_ray_xyz_with_gradient(ray_in_xyz, face_norm, n);

for i = 1:num
    curr_rot_mat = rotz(90 + rot_llq(1)) * rotx(90 - rot_llq(2)) * rotz(rot_llq(3));
    ray_out_xyz(i, :) = ray_out_xyz(i, :) * curr_rot_mat';
end
ray_out_ll = xyz2ll_with_gradient(ray_out_xyz);
end