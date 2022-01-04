function [ray_out_ll, jac] = crystal_system(rot, ray_in_ll, crystal, trace)
% Find output ray of a crystal system with given input ray and crystal rotation.
%
% INPUT
%   rot:                n*3 or n*4, it has either of the tow form (can be specified explicitly):
%                       1. [longitude, latitude, roll], in degree. A representation of crystal rotation.
%                          see help of geo.llr2quat(llr) for detail. This form is default.
%                       2. [w, x, y, z] for quaternion [w, xi, yj, zk].
%   ray_in_ll:          1*2, [longitude, latitude], in degree. Input ray direction.
%   crystal:            struct
%   trace:              struct
%
% OUTPUT
%   ray_out_ll:         n*2, [longitude, latitude], in degree. Output ray direction.
%   jac:                2*3*n, Jacobian. Input is rotation (llr or quat3) and output is ray_out_ll

num = size(rot, 1);
dim_rot = size(rot, 2);

xyz0 = geo.ll2xyz(ray_in_ll);
if dim_rot == 3
    [rot_mat, jac_rot_mat] = geo.llr2mat(rot);
elseif dim_rot == 4
    [rot_mat, jac_rot_mat] = geo.quat2mat(rot);
else
    error('rot must have 3 or 4 columns!');
end

ray_in_xyz = zeros(num, 3);
for i = 1:num
    ray_in_xyz(i, :) = xyz0 * rot_mat(:, :, i)';
end

[ray_out_xyz, jac_out_xyz] = opt.trace_ray_xyz(ray_in_xyz, crystal, trace);

jac_xyz = zeros(3, dim_rot, num);
for i = 1:num
    tmp_xyz = ray_out_xyz(i, :)';
    if dim_rot == 3
        g_xyz1 = [jac_rot_mat(:, :, 1, i) * xyz0', jac_rot_mat(:, :, 2, i) * xyz0', ...
                jac_rot_mat(:, :, 3, i) * xyz0'];
        jac_xyz(:, :, i) = [jac_rot_mat(:, :, 1, i)' * tmp_xyz, jac_rot_mat(:, :, 2, i)' * tmp_xyz, ...
                            jac_rot_mat(:, :, 3, i)' * tmp_xyz] + ...
            rot_mat(:, :, i)' * jac_out_xyz(:, :, i) * g_xyz1;
    else
        g_xyz1 = [jac_rot_mat(:, :, 1, i) * xyz0', jac_rot_mat(:, :, 2, i) * xyz0', ...
                jac_rot_mat(:, :, 3, i) * xyz0', jac_rot_mat(:, :, 4, i) * xyz0'];
        jac_xyz(:, :, i) = [jac_rot_mat(:, :, 1, i)' * tmp_xyz, jac_rot_mat(:, :, 2, i)' * tmp_xyz, ...
                            jac_rot_mat(:, :, 3, i)' * tmp_xyz, jac_rot_mat(:, :, 4, i)' * tmp_xyz] + ...
            rot_mat(:, :, i)' * jac_out_xyz(:, :, i) * g_xyz1;
    end
    ray_out_xyz(i, :) = ray_out_xyz(i, :) * rot_mat(:, :, i);
end

[ray_out_ll, jac_out_ll] = geo.xyz2ll(ray_out_xyz);
jac = zeros(2, dim_rot, num);
for i = 1:num
    jac(:, :, i) = jac_out_ll(:, :, i) * jac_xyz(:, :, i);
end
end
