function config = init_config_3d(crystal, trace, sun_ll, level)
% initial grid
[~, rot0_ll, dr] = generate_healpix_grids(level);
roll0 = (0:dr*2:360)';

a_num = size(rot0_ll, 1);
r_num = length(roll0);
total_num = a_num * r_num;

axis_llr_store = [repmat(rot0_ll, [r_num, 1]), kron(roll0, ones(a_num, 1))];

out_xyz = zeros(total_num, 3);
out_ll = zeros(total_num, 2);
mat_store = zeros(2, 3, total_num);

progress_cnt = 0;
progress_bin = 0.1;
for i = 1:total_num
    if progress_cnt > progress_bin
        fprintf('init_config_3d (%04.1f%%)...\n', i / total_num * 100);
        progress_cnt = progress_cnt - floor(progress_cnt / progress_bin) * progress_bin;
    end
    progress_cnt = progress_cnt + 1 / total_num;

    [tmp_out_ll, tmp_jacob] = opt.crystal_system(axis_llr_store(i, :), sun_ll, crystal, trace);
    out_ll(i, :) = tmp_out_ll;
    out_xyz(i, :) = geo.ll2xyz(tmp_out_ll);
    mat_store(:, :, i) = tmp_jacob;
end

sun_xyz = geo.ll2xyz(sun_ll);
bending_angle = acosd(out_xyz * sun_xyz');

config.bending_angle = bending_angle;
config.axis_llr_store = axis_llr_store;
config.axis_quat_store = geo.llr2quat(axis_llr_store);
config.dr = dr;
config.out_ll = out_ll;
config.out_xyz = out_xyz;
config.jacobian = mat_store;
end