function test_find_contour_fdf()
% Test function find_contour_fdf

test_cases = {@suite1};
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
for i = 1:size(x, 1)
    assert(abs(fdf1(x(i, :)) - y0) < 1e-6);
end
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
