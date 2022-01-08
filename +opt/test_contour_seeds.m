function test_contour_seeds()
% Test a real case for contour seeds.

test_cases = {@suite1};
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

% Start
[~, seed_quat, ~] = opt.find_seed_rot(config, ray_out_ll);

seed_rot = seed_quat;
contour_h = 0.05;
reduce_eps = 0.03;

fprintf(' case 1 ... ');
rot0 = [-0.516278096307973, 0.170107948397147, -0.254138559784028, -0.799958627442999];
[rot_contour, ~] = ode.find_contour_fdf(fdf, rot0, 'h', contour_h);
[seed_rot, ~] = geo.reduce_pts_polyline(rot_contour, seed_rot, 'eps', reduce_eps);
if debug
    figure(1); clf;
    hold on;
    plot_rot_space(seed_quat, [2, 3, 4], 'o');
    plot_rot_space(seed_rot, [2, 3, 4], 's');
    plot_rot_space(rot_contour, [2, 3, 4], '-o');
    plot_rot_tan_space(fdf, rot_contour, [2, 3, 4], 'm');
    axis equal;
end
fprintf('passed!\n');

fprintf(' case 2 ... ');
rot0 = [0.588369809379402, -0.347408255513868, 0.683166656311173, -0.257704852721616];
[rot_contour, ~] = ode.find_contour_fdf(fdf, rot0, 'h', contour_h);
[seed_rot, ~] = geo.reduce_pts_polyline(rot_contour, seed_rot, 'eps', reduce_eps);
if debug
    figure(2); clf;
    hold on;
    plot_rot_space(seed_quat, [2, 3, 4], 'o');
    plot_rot_space(seed_rot, [2, 3, 4], 's');
    plot_rot_space(rot_contour, [2, 3, 4], '-o');
    plot_rot_tan_space(fdf, rot_contour, [2, 3, 4], 'm');
    axis equal;
end
fprintf('passed!\n');
end

% ================================================================================
function plot_rot_space(rot, space_idx, varargin)
if isempty(space_idx)
    space_idx = [1, 2, 3];
end
plot3(rot(:, space_idx(1)), rot(:, space_idx(2)), rot(:, space_idx(3)), varargin{:});
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
