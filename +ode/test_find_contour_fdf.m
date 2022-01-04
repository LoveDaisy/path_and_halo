function test_find_contour_fdf()
% Test function find_contour_fdf

test_cases = {@suite1, @suite2, @suite3};
num = length(test_cases);

fprintf('Start testing for find_contour_fdf...\n');
for i = 1:num
    fprintf('Testing suite %d/%d...\n', i, num);
    test_cases{i}();
    fprintf('suite %d/%d passed!\n', i, num);
end
end

% ================================================================================
function suite1()
x0 = [0.5, 0.5];
y0 = fdf1(x0);
[x, status] = ode.find_contour_fdf(@fdf1, x0);

assert(status.closed);
assert(status.completed);
assert(all(abs(fdf1(x) - y0) < 1e-6));
end

% ================================================================================
function suite2()
fprintf(' case 1 ... ');
x0 = [0.6, -1];
y0 = fdf2(x0);
[x, status] = ode.find_contour_fdf(@fdf2, x0);

assert(~status.closed);
assert(status.completed);
assert(all(abs(fdf2(x) - y0) < 1e-6));
assert(abs(min(x(:, 1))) < 0.01);
assert(abs(max(x(:, 1)) - 1) < 0.01);
fprintf('passed!\n');

fprintf(' case 2 ... ');
x0 = [0.5, 0];
y0 = fdf2(x0);
[x, status] = ode.find_contour_fdf(@fdf2, x0);

assert(~status.closed);
assert(status.completed);
assert(all(abs(fdf2(x) - y0) < 1e-6));
assert(abs(min(x(:, 1))) < 0.01);
assert(abs(max(x(:, 1)) - 1) < 0.01);
fprintf('passed!\n');
end

% ================================================================================
function suite3()
crystal = opt.make_prism_crystal(1);
trace.fid = [3; 5];
sun_ll = [180, 25];
ray_in_ll = [sun_ll(1) + 180, -sun_ll(2)];

config = opt.init_config(crystal, trace, sun_ll, 3);

fprintf(' case 1 ... ');
fdf = @(rot) opt.crystal_system(rot, ray_in_ll, crystal, trace, 'llr');

ray_out_ll = [0, -1];
ray_out_xyz = geo.ll2xyz(ray_out_ll);

init_angle_diff = acosd(config.out_xyz * ray_out_xyz');
[~, idx] = min(init_angle_diff);
init_llr = ode.find_solution_fdf(fdf, config.axis_llr_store(idx, :), ray_out_ll);
[rot_llr, status] = ode.find_contour_fdf(fdf, init_llr, 'eps', 1e-4, 'h', 0.5);
assert(status.closed);
for i = 1:size(rot_llr, 1)
    assert(norm(fdf(rot_llr(i, :)) - ray_out_ll) < 1e-4);
end
fprintf('passed!\n');
end

% ================================================================================
function [y, jac] = fdf1(x)
num = size(x, 1);

y = exp(x(:, 1) + 3 * x(:, 2) - 0.1) + exp(x(:, 1) - 3 * x(:, 2) - 0.1) + exp(-x(:, 1) - 0.1);
jac = zeros(1, 2, num);
for i = 1:num
    jac(:, :, i) = [exp(x(i, 1) + 3 * x(i, 2) - 0.1) + exp(x(i, 1) - 3 * x(i, 2) - 0.1) - exp(-x(i, 1) - 0.1), ...
                    3 * exp(x(i, 1) + 3 * x(i, 2) - 0.1) - 3 * exp(x(i, 1) - 3 * x(i, 2) - 0.1)];
end
end

% ================================================================================
function [y, jac] = fdf2(x)
num = size(x, 1);

y = log((x(:, 1).^2 - x(:, 2)).^2) - log(x(:, 1)) - log(1 - x(:, 1));
jac = [4 * x(:, 1) ./ (x(:, 1).^2 - x(:, 2)) + 1 ./ (1 - x(:, 1)) - 1 ./ x(:, 1), ...
        - 2 ./ (x(:, 1).^2 - x(:, 2))]';
jac = reshape(jac, [1, 2, num]);
end
