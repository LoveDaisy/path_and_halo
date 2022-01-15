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
p.addParameter('MaxEval', 30, @(x) validateattributes(x, {'double'}, {'scalar', 'integer', 'positive'}));
p.parse(varargin{:});

[y0, jac] = fun(x0);
h_min = 0.01;

dy = yq - y0;
x = x0;
fun_eval_cnt = 1;
h = 1;
while norm(dy) > p.Results.eps && fun_eval_cnt < p.Results.MaxEval && h > h_min
    % Find deepest gradient
    dx = jac \ dy';
    jn = null(jac);
    dx = dx - dot(dx, jn(:, 1)) * jn;

    % Linear search
    h = 1;
    x = x0 + (dx .* h)';
    [y, jac1] = fun(x);
    fun_eval_cnt = fun_eval_cnt + 1;

    convexity = eps_sign((y - y0)' - jac * (dx .* h), 1e-4);
    k_direction = eps_sign(jac * dx, 1e-4);

    a = 4.^(convexity .* k_direction);
    b = 0.6;
    indicating_flag = (convexity .* (abs(k_direction) > 0.5));
    linearity = (y - y0)' - a .* (jac * (dx .* h));
    h_shrink_idx = eps_sign(linearity .* indicating_flag, 1e-4) > 0.5;
    while norm(yq - y) > p.Results.eps && fun_eval_cnt < p.Results.MaxEval && ...
            h > h_min && any(h_shrink_idx)
        h = h * b;
        x = x0 + (dx .* h)';
        [y, jac1] = fun(x);
        fun_eval_cnt = fun_eval_cnt + 1;

        convexity = eps_sign((y - y0)' - jac * (dx .* h), 1e-4);
        k_direction = eps_sign(jac * dx, 1e-4);
        a = 4.^(convexity .* k_direction);
        indicating_flag = (convexity .* (abs(k_direction) > 0.5));
        linearity = (y - y0)' - a .* (jac * (dx .* h));
        h_shrink_idx = eps_sign(linearity .* indicating_flag, 1e-4) > 0.5;
    end

    x0 = x;
    y0 = y;
    dy = yq - y;
    jac = jac1;
end
status.finish = norm(dy) < p.Results.eps;
status.fun_eval_cnt = fun_eval_cnt;
end

function y = eps_sign(x, e)
y = sign(x);
y(abs(x) < e) = 0;
end
