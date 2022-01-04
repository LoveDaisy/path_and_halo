function [x, status] = find_solution_fdf(fdf, x0, yq, varargin)
% Find the exact solution for a fdf function.
% An fdf function is a function that gives both value and gradient (Jacobian)
%
% INPUT
%  fdf:         function that outputs both value (y) and gradient (or Jacobian, jac). [y, jac] = fun(x)
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
p.addParameter('MaxEval', 10, @(x) validateattributes(x, {'double'}, {'scalar', 'integer', 'positive'}));
p.parse(varargin{:});

fun = @(x) jac_warpper(fdf, x);
[y0, jac] = fun(x0);

dy = yq - y0;
x = x0;
fun_eval_cnt = 1;
while norm(dy) > p.Results.eps && fun_eval_cnt < p.Results.MaxEval
    % Find deepest gradient
    dx = jac \ dy';

    % Linear search
    h = 1;
    x = x0 + dx' * h;
    [y, jac1] = fun(x);
    fun_eval_cnt = fun_eval_cnt + 1;

    convexity = ((y - y0)' > h * jac * dx) * 2 - 1;
    k_direction = (jac * dx > 0) * 2 - 1;

    a = 4.^(convexity .* k_direction);
    b = 0.6;
    while norm(yq - y) > p.Results.eps && fun_eval_cnt < p.Results.MaxEval && ...
            h > 0.01 && any((y - y0)' .* convexity > (h * a .* convexity) .* (jac * dx))
        h = h * b;
        x = x0 + dx' * h;
        [y, jac1] = fun(x);
        fun_eval_cnt = fun_eval_cnt + 1;

        convexity = ((y - y0)' > h * jac * dx) * 2 - 1;
        a = 4.^(convexity .* k_direction);
    end

    x0 = x;
    y0 = y;
    dy = yq - y;
    jac = jac1;
end
status.finish = norm(dy) < p.Results.eps;
status.fun_eval_cnt = fun_eval_cnt;
end

function [y, jac] = jac_warpper(fun, x)
[y, jac] = fun(x);
dim2 = size(y, 2);
if size(jac, 1) > dim2
    jac = jac(1:dim2, :);
end
end
