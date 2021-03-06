function test_contour_seeds()
% Test a real case for contour seeds.

test_cases = {@suite1, @suite2};
num = length(test_cases);
debug = true;

fprintf('Start testing for contour seeds...\n');
for i = 1:num
    fprintf('Testing suite %d/%d...\n', i, num);
    test_cases{i}(debug);
    fprintf('suite %d/%d passed!\n', i, num);
end
close all;
end

% ================================================================================
function suite1(debug)
crystal = opt.make_prism_crystal(1);
trace.fid = [3; 5];
sun_ll = [180, 10];
ray_in_ll = [sun_ll(1) + 180, -sun_ll(2)];
ray_out_ll = [-5, 13.8];
ray_out_xyz = geo.ll2xyz(ray_out_ll);

config_cache_file = 'test_config_1_35_180+10_3.mat';
if exist(config_cache_file, 'file')
    load(config_cache_file);
else
    config = opt.init_config(crystal, trace, sun_ll);
    save(config_cache_file, 'config');
end

fdf = @(rot) opt.crystal_system(rot, ray_in_ll, crystal, trace);

% Find seed rotation
cand_rot = opt.find_cand_rot(config, ray_out_ll, 'quat');
init_rot = ode.find_solution(fdf, cand_rot(1, :), [ray_out_xyz, 1], 'eps', 1e-8);

contour_h = 0.05;
reduce_eps = config.dr * 1.5;

% Find contour
rot_contour = ode.find_contour(fdf, init_rot, 'h', contour_h);

% Reduce seed points
seed_rot = geo.reduce_pts_polyline(rot_contour, cand_rot, 'eps', reduce_eps, 'jac_fun', fdf);

if debug
    space_idx = [2, 3, 4];
    figure(1); clf;
    hold on;
    utl.plot_data_3d(cand_rot, space_idx, 'o');
    utl.plot_data_3d(seed_rot, space_idx, 's');
    utl.plot_data_3d(rot_contour, space_idx, '-x');
    utl.plot_tan_space_3d(rot_contour, fdf, space_idx);
    axis equal;
end

assert(size(seed_rot, 1) / size(cand_rot, 1) < 0.1);
end

% ================================================================================
function suite2(debug)
crystal = opt.make_prism_crystal(1);
trace.fid = [1; 3; 2; 4; 5; 1];
sun_ll = [0, 10];
ray_in_ll = [sun_ll(1) + 180, -sun_ll(2)];

config_cache_file = 'test_config_1_132451_0-10_3.mat';
if exist(config_cache_file, 'file')
    load(config_cache_file);
else
    config = opt.init_config(crystal, trace, sun_ll);
    save(config_cache_file, 'config');
end

contour_h = 0.05;
reduce_eps = 0.03;
fdf = @(rot) opt.crystal_system(rot, ray_in_ll, crystal, trace);

ray_out_ll = [0, 10];
cand_rot = opt.find_cand_rot(config, ray_out_ll, 'quat');
seed_rot = cand_rot;

% --------
fprintf(' case 1 ... ');
rot1 = [-0.625137188553982, -0.0773308939540198, 0.607468716849977, 0.483947503748253];
[rot_contour1, ~] = ode.find_contour(fdf, rot1, 'h', contour_h);
[seed_rot, ~] = geo.reduce_pts_polyline(rot_contour1, seed_rot, 'eps', reduce_eps);
if debug
    space_idx = [2, 3, 4];
    figure(1); clf;
    hold on;
    utl.plot_data_3d(cand_rot, space_idx, 'o');
    utl.plot_data_3d(seed_rot, space_idx, 's');
    utl.plot_data_3d(rot_contour1, space_idx, '-x');
    utl.plot_tan_space_3d(rot_contour1, fdf, space_idx);
    axis equal;
end
fprintf('passed!\n');

% --------
fprintf(' case 2 ... ');
rot1 = [-0.203879731194181, -0.489429488302747, 0.847848784604723, 0.00665354525130694];
rot_contour1 = ode.find_contour(fdf, rot1, 'h', contour_h);
seed_rot = geo.reduce_pts_polyline(rot_contour1, seed_rot, 'eps', reduce_eps);
if debug
    space_idx = [2, 3, 4];
    figure(2); clf;
    hold on;
    utl.plot_data_3d(cand_rot, space_idx, 'o');
    utl.plot_data_3d(seed_rot, space_idx, 's');
    utl.plot_data_3d(rot_contour1, space_idx, '-x');
    utl.plot_tan_space_3d(rot_contour1, fdf, space_idx);
    axis equal;
end
fprintf('passed!\n');

% --------
ray_out_ll = [0.5, 10];
ray_out_xyz = geo.ll2xyz(ray_out_ll);
cand_rot = opt.find_cand_rot(config, ray_out_ll, 'quat');
init_rot = ode.find_solution(fdf, cand_rot(1, :), [ray_out_xyz, 1], 'eps', 1e-8);

fprintf(' case 3 ... ');
rot_contour1 = ode.find_contour(fdf, init_rot, 'h', contour_h);
seed_rot = geo.reduce_pts_polyline(rot_contour1, seed_rot, 'eps', reduce_eps);
if debug
    space_idx = [2, 3, 4];
    figure(3); clf;
    hold on;
    utl.plot_data_3d(cand_rot, space_idx, 'o');
    utl.plot_data_3d(seed_rot, space_idx, 's');
    utl.plot_data_3d(rot_contour1, space_idx, '-x');
    utl.plot_tan_space_3d(rot_contour1, fdf, space_idx);
    axis equal;
end
fprintf('passed!\n');

% --------
ray_out_ll = [0.5, 7.5];
ray_out_xyz = geo.ll2xyz(ray_out_ll);
seed_quat = opt.find_cand_rot(config, ray_out_ll, 'quat');

fprintf(' case 4 ... ');
rot1 = ode.find_solution(fdf, seed_quat(1, :), [ray_out_xyz, 1], 'eps', 1e-8);
rot_contour1 = ode.find_contour(fdf, rot1, 'h', contour_h);
seed_quat1 = geo.reduce_pts_polyline(rot_contour1, seed_quat, 'jac_fun', fdf, 'eps', config.dr);
if debug
    space_idx = [2, 3, 4];
    figure(4); clf;
    hold on;
    utl.plot_data_3d(seed_quat, space_idx, 'o');
    utl.plot_data_3d(seed_quat1, space_idx, 'x');
    utl.plot_data_3d(rot_contour1, space_idx, '-^');
    utl.plot_data_3d(seed_quat(1, :), space_idx, 's', 'markersize', 15);
    utl.plot_tan_space_3d(rot_contour1, fdf, space_idx);
    axis equal;
end

rot2 = ode.find_solution(fdf, seed_quat(end, :), [ray_out_xyz, 1], 'eps', 1e-8);
rot_contour2 = ode.find_contour(fdf, rot2, 'h', contour_h);
seed_quat2 = geo.reduce_pts_polyline(rot_contour2, seed_quat, 'jac_fun', fdf, 'eps', config.dr);
if debug
    space_idx = [2, 3, 4];
    figure(5); clf;
    hold on;
    utl.plot_data_3d(seed_quat, space_idx, 'o');
    utl.plot_data_3d(seed_quat2, space_idx, 'x');
    utl.plot_data_3d(rot_contour2, space_idx, '-^');
    utl.plot_data_3d(seed_quat(end, :), space_idx, 's', 'markersize', 15);
    utl.plot_tan_space_3d(rot_contour2, fdf, space_idx);
    axis equal;
end
fprintf('passed!\n');
end
