function test_find_solution_fdf()
test_fdf = {@case1, @case2};
test_x0 = {[1.2, 1.1], [0.5, 0.5]};
test_yq = {2, 4.5};

test_cases = struct('fdf', test_fdf, 'x0', test_x0, 'yq', test_yq);
case_num = length(test_cases);

fprintf('Start tesing find_solution_fdf...\n');
tol_eps = 1e-8;
for i = 1:case_num
    fprintf('Testing case %d/%d...', i, case_num);
    x0 = test_cases(i).x0;
    yq = test_cases(i).yq;
    fdf = test_cases(i).fdf;
    [x, ~, num] = ode.find_solution_fdf(fdf, x0, yq, 'eps', tol_eps);
    y = fdf(x);
    assert(norm(y - yq) < tol_eps, 'Test case %d failed!', i);
    fprintf(' passed! Function calls: %d\n', num);
end
end

function [y, jac] = case1(x)
%  y = x1^2 + x2^2

num = size(x, 1);

y = sum(x.^2, 2);
jac = 2 * reshape(x', [1, 2, num]);
end

function [y, jac] = case2(x)
%  y = exp(x1 + 3*x2 - 0.1) + exp(x1 - 3*x2 - 0.1) + exp(-x1 - 0.1);
num = size(x, 1);

y = exp(x(:, 1) + 3 * x(:, 2) - 0.1) + exp(x(:, 1) - 3 * x(:, 2) - 0.1) + exp(-x(:, 1) - 0.1);
jac = zeros(1, 2, num);
for i = 1:num
    jac(:, :, i) = [exp(x(i, 1) + 3 * x(i, 2) - 0.1) + exp(x(i, 1) - 3 * x(i, 2) - 0.1) - exp(-x(i, 1) - 0.1), ...
                    3 * exp(x(i, 1) + 3 * x(i, 2) - 0.1) - 3 * exp(x(i, 1) - 3 * x(i, 2) - 0.1)];
end
end
