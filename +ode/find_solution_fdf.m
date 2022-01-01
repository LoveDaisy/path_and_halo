function [x, flag, fun_eval_num] = find_solution_fdf(fun, x0, yq, varargin)
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
%  flag:            bool, whether find a solution within given tolerance.
%  fun_eval_num:    how many times the function is called.

num = size(x0, 1);
dim1 = size(x0, 2);
dim2 = size(yq, 2);

p = inputParser;
p.addRequired('fun', @(x) validateattributes(x, {'function_handle'}, {'scalar'}));
p.addRequired('x0', @(x) validateattributes(x, {'numeric'}, {'2d'}));
p.addRequired('yq', @(x) validateattributes(x, {'numeric'}, {'2d', 'nrows', num}));
p.addParameter('eps', 1e-8, @(x) validateattributes(x, {'double'}, {'scalar', 'positive'}));
p.addParameter('MaxEval', 10, @(x) validateattributes(x, {'double'}, {'scalar', 'integer', 'positive'}));
p.parse(fun, x0, yq, varargin{:});

[y0, jac] = fun(x0);
validateattributes(y0, {'numeric'}, {'2d', 'size', [num, dim2]});
validateattributes(jac, {'numeric'}, {'3d', 'size', [dim2, dim1, num]});

dy = yq - y0;
fun_eval_num = 1;
while norm(dy) > p.Results.eps && fun_eval_num < p.Results.MaxEval
    % Find deepest gradient
    dx = jac \ dy;

    % Linear search
    h = 1;
    x = x0 + dx' * h;
    [y, jac] = fun(x);
    fun_eval_num = fun_eval_num + 1;

    a = 0.25;
    b = 0.6;
    while norm(yq - y) > p.Results.eps && fun_eval_num < p.Results.MaxEval && ...
            h > 0.01 && any(y - y0 > a * h * jac * dx)
        h = h * b;
        x = x0 + dx' * h;
        [y, jac] = fun(x);
        fun_eval_num = fun_eval_num + 1;
    end

    x0 = x;
    y0 = y;
    dy = yq - y;
end
flag = norm(dy) < p.Results.eps;
end
