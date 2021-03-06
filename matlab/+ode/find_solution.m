function [x, status] = find_solution(fun, x0, yq, varargin)
% Find the exact solution for a fdf function.
% An fdf function is a function that gives both value and gradient (Jacobian)
%
% INPUT
%  fun:         function that outputs both value (y) and gradient (or Jacobian, jac). [y, jac] = fun(x)
%  x0:          n*d1, start point, fun(x0) should be approximate to yq
%  yq:          n*d2, target value
%
% OPTIONS
%  'eps':       absolute tolerance for output. norm(fun(x) - yq) < eps. default 1e-8
%  'MaxEval':   maximumn number for function evaluation
%
% OUTPUT
%  x:               n*d1, exact solution
%  status:          struct, contains some information as following:
%    .finish:       bool, whether find a solution within given tolerance.
%    .fun_eval_cnt: how many times the function is called.

p = inputParser;
p.addParameter('eps', 1e-8, @(x) validateattributes(x, {'double'}, {'scalar', 'positive'}));
p.addParameter('MaxEval', 15, @(x) validateattributes(x, {'double'}, {'scalar', 'integer', 'positive'}));
p.parse(varargin{:});

status.finish = false;
status.fun_eval_cnt = 0;

[y0, jac] = fun(x0);
if any(isnan(y0))
    x = x0;
    return;
end
h_min = 0.1;

dy = yq - y0;
x = x0;
fun_eval_cnt = 1;
h = 1;
while norm(dy) > p.Results.eps && fun_eval_cnt < p.Results.MaxEval && h > h_min
    % Find deepest gradient
    r_jac = rank(jac, 1e-8);
    [u, s, v] = svd(jac);
    dx = v(:, 1:r_jac) / s(1:r_jac, 1:r_jac) * u(:, 1:r_jac)' * dy';

    % Linear search
    h = 1;
    x = x0 + (dx .* h)';
    y = fun(x);
    fun_eval_cnt = fun_eval_cnt + 1;

    val_nan = any(isnan(y));
    b = 0.6;
    while any(isnan(y)) || (norm(yq - y) > p.Results.eps && fun_eval_cnt < p.Results.MaxEval && ...
            h > h_min && norm(y - yq) > norm(dy) * h * 0.8)
        h = h * b;
        x = x0 + (dx .* h)';
        y = fun(x);
        fun_eval_cnt = fun_eval_cnt + 1;
    end

    if (val_nan && h < 0.3) || (h < h_min)
        break;
    end

    x0 = x;
    dy = yq - y;
    [~, jac] = fun(x);
end
status.finish = norm(dy) < p.Results.eps;
status.fun_eval_cnt = fun_eval_cnt;
end

function y = eps_sign(x, e)
y = sign(x);
y(abs(x) < e) = 0;
end
