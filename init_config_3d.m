function config = init_config_3d(face_norm, n, sun_ll, level)
n = [1; n(:)];

% initial grid
[~, rot0_ll, dr] = generate_healpix_grids(level);
roll0 = (0:dr*2:360)';

a_num = size(rot0_ll, 1);
r_num = length(roll0);
total_num = a_num * r_num;

axis_rot_store = [repmat(rot0_ll, [r_num, 1]), kron(roll0, ones(a_num, 1))];

out_xyz = zeros(total_num, 3);
out_ll = zeros(total_num, 2);
mat_store = zeros(3, 3, total_num);
for i = 1:total_num
    [tmp_out_ll, tmp_jacob] = crystal_system_with_gradient(axis_rot_store(i, :), sun_ll, face_norm, n);
    out_ll(i, :) = tmp_out_ll;
    out_xyz(i, :) = ll2xyz_with_gradient(tmp_out_ll);
    mat_store(:, :, i) = tmp_jacob;
end

sun_xyz = ll2xyz_with_gradient(sun_ll);
bending_angle = acosd(sum(sun_xyz .* out_xyz, 2));
bending_angle_max = max(bending_angle);
bending_angle_min = min(bending_angle);

config.bending_angle = bending_angle;
config.bending_angle_max = bending_angle_max;
config.bending_angle_min = bending_angle_min;
config.axis_rot_store = axis_rot_store;
config.rot0_ll = rot0_ll;
config.roll0 = roll0;
config.dr = dr;
config.out_ll = out_ll;
config.out_xyz = out_xyz;
config.mat_store = mat_store;
end