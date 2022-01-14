function test_contour_seeds()
% Test a real case for contour seeds.

test_cases = {@suite1, @suite2};
num = length(test_cases);
debug = false;

fprintf('Start testing for contour seeds...\n');
for i = 1:num
    fprintf('Testing suite %d/%d...\n', i, num);
    test_cases{i}(debug);
    fprintf('suite %d/%d passed!\n', i, num);
end
end

% ================================================================================
function suite1(debug)
crystal = opt.make_prism_crystal(1);
trace.fid = [3; 5];
sun_ll = [180, 10];
ray_in_ll = [sun_ll(1) + 180, -sun_ll(2)];
ray_out_ll = [-5, 13.8];

config_cache_file = 'test_config_1_35_180+10_3.mat';
if exist(config_cache_file, 'file')
    load(config_cache_file);
else
    config = opt.init_config(crystal, trace, sun_ll, 3);
    save(config_cache_file, 'config');
end

fdf = @(rot) opt.crystal_system(rot, ray_in_ll, crystal, trace);

% Find seed rotation
[seed_quat, ~] = opt.find_seed_rot(config, ray_out_ll, 'quat');

seed_rot = seed_quat;
contour_h = 0.05;
reduce_eps = 0.03;

% Find contour
[rot_contour, ~] = ode.find_contour(fdf, seed_rot(1, :), 'h', contour_h);

% Reduce seed points
[seed_rot, ~] = geo.reduce_pts_polyline(rot_contour, seed_rot, 'eps', reduce_eps);

if debug
    space_idx = [2, 3, 4];
    figure(1); clf;
    hold on;
    plot_data_3d(seed_quat, space_idx, 'o');
    plot_data_3d(seed_rot, space_idx, 's');
    plot_data_3d(rot_contour, space_idx, '-o');
    plot_rot_tan_space(fdf, rot_contour, space_idx, 'm');
    axis equal;
end

assert(isempty(seed_rot));
end

% ================================================================================
function suite2(debug)
crystal = opt.make_prism_crystal(1);
trace.fid = [1; 3; 2; 4; 5; 1];
sun_ll = [0, 10];
ray_in_ll = [sun_ll(1) + 180, -sun_ll(2)];
ray_out_ll = [0, 10];

config_cache_file = 'test_config_1_132451_0-10_3.mat';
if exist(config_cache_file, 'file')
    load(config_cache_file);
else
    config = opt.init_config(crystal, trace, sun_ll, 3);
    save(config_cache_file, 'config');
end

fdf = @(rot) opt.crystal_system(rot, ray_in_ll, crystal, trace);

[seed_quat, ~] = opt.find_seed_rot(config, ray_out_ll, 'quat');

seed_rot = seed_quat;
contour_h = 0.05;
reduce_eps = 0.03;

fprintf(' case 1 ... ');
rot0 = [-0.625137188553982,-0.0773308939540198,0.607468716849977,0.483947503748253];
[rot_contour, ~] = ode.find_contour(fdf, rot0, 'h', contour_h);
[seed_rot, ~] = geo.reduce_pts_polyline(rot_contour, seed_rot, 'eps', reduce_eps);
if debug
    space_idx = [2, 3, 4];
    figure(1); clf;
    hold on;
    plot_data_3d(seed_quat, space_idx, 'o');
    plot_data_3d(seed_rot, space_idx, 's');
    plot_data_3d(rot_contour, space_idx, '-o');
    plot_rot_tan_space(fdf, rot_contour, space_idx, 'm');
    axis equal;
end
fprintf('passed!\n');

fprintf(' case 2 ... ');
rot0 = [-0.203879731194181,-0.489429488302747,0.847848784604723,0.00665354525130694];
[rot_contour, ~] = ode.find_contour(fdf, rot0, 'h', contour_h);
[seed_rot, ~] = geo.reduce_pts_polyline(rot_contour, seed_rot, 'eps', reduce_eps);
if debug
    space_idx = [2, 3, 4];
    figure(2); clf;
    hold on;
    plot_data_3d(seed_quat, space_idx, 'o');
    plot_data_3d(seed_rot, space_idx, 's');
    plot_data_3d(rot_contour, space_idx, '-o');
    plot_rot_tan_space(fdf, rot_contour, space_idx, 'm');
    axis equal;
end
fprintf('passed!\n');
end

% ================================================================================
function plot_rot_tan_space(fdf, rot, space_idx, varargin)
if isempty(space_idx)
    space_idx = [1, 2, 3];
end

num = size(rot, 1);

h = 5;
hold on;
for i = 1:num
    [~, jac] = fdf(rot(i, :));

    jjt = jac(1:2, :) * jac(1:2, :)';
    jjt2 = jjt * jjt;

    cos_q = cosd(0:10:360)';
    sin_q = sind(0:10:360)';

    v2 = sqrt(jjt2(1, 1) / max(jjt2(1, 1) * jjt2(2, 2) - jjt2(1, 2)^2, 1e-3)) * sin_q;
    v1 = cos_q / sqrt(jjt2(1, 1)) - jjt2(1, 2) / jjt2(1, 1) * v2;
    u = [v1, v2] * jac(1:2, :) * h;
    u = bsxfun(@plus, u, rot(i, :));
    plot3(u(:, space_idx(1)), u(:, space_idx(2)), u(:, space_idx(3)), varargin{:});
end
end
