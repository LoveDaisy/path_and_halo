function test_compute_weight()
test_cases = {@suite1};
num = length(test_cases);
debug = false;

fprintf('Start testing for compute_weight...\n');
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

% Find contour
[rot_contour, contour_status] = ode.find_contour(fdf, seed_quat(1, :), 'h', 0.05);
assert(contour_status.closed && contour_status.completed);

rot_llr = geo.quat2llr(rot_contour);
w0 = 0.00569990578729381;

axis_pdf = generate_axis_pdf([0, 0, 0]);
[w, cmp] = opt.compute_contour_weight(rot_llr, axis_pdf, config);

assert(abs(w - w0) < 1e-8);
end

% ================================================================================
function axis_pdf = generate_axis_pdf(zenith_dist)
% INPUT
%   zenith_dist:        [type, mean, sigma], type: 0 - uniform, 1 - gaussian

if zenith_dist(1) == 1
    zenith_dist = zenith_dist(2:3);
    tmp_x = linspace(-90, 90, 50000);
    tmp_pdf = exp(- (90 - tmp_x - zenith_dist(1)).^2/2 / zenith_dist(2)^2) / zenith_dist(2);
    zen_total = (sum(tmp_pdf) * (tmp_x(2) - tmp_x(1)));
    axis_pdf = @(llr) (exp(- (90 - asind(sind(llr(:, 2))) - zenith_dist(1)).^2/2 / zenith_dist(2)^2) / ...
        zenith_dist(2) / zen_total) * (1/360) * (1/360);
elseif zenith_dist(1) == 0
    axis_pdf = @(llr) ones(size(llr, 1), 1) / (180 * 180 * 4 / pi);
else
    error('Distribution type error!');
end
end
