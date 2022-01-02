function [x, status] = find_contour_fdf(fdf, x0, varargin)
% Find a contour for an ODE with given target value.
% For every point xi of output x, fdf(xi) == fdf(x0).
%
% INPUT
%   fdf:        ODE function, [y, jac] = fdf(x), input x and output y may both be vector
%   x0:         scalar. Initial point.
%
% OPTION
%   'eps':      scalar. Tolerance for output y.
%   'MaxPts':   maximumn points in the result contour. Default 100
%
% OUTPUT
%   x:              n*d
%   status:         struct, state when exit. Including following fields:
%     .closed       bool, whether this contour is a closed loop
%     .completed    bool, whether this contour is completed, i.e. both ends cannot be extended
%     .fun_eval_cnt integer, how many times that function is called

p = inputParser;
p.addRequired('fdf', @(x) validateattributes(x, {'function_handle'}, {'scalar'}));
p.addRequired('x0', @(x) validateattributes(x, {'numeric'}, {'nrows', 1}));
p.addParameter('eps', 1e-8, @(x) validateattributes(x, {'double'}, {'scalar', 'positive'}));
p.addParameter('MaxPts', 100, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
p.parse(fdf, x0, varargin{:});

h = 0.05;
search_options = {'method', 'euler', 'eps', p.Results.eps, 'MaxPts', p.Results.MaxPts, 'h', h};
% Search forward direction
[x_f, status_f] = search_direction(fdf, x0, 1, search_options{:});
status = status_f;
if status_f.closed
    x = xf;
    return;
end

% Search backward direction
[x_b, status_b] = search_direction(fdf, x0, -1, search_options{:});
if isempty(x_b)
    return;
end

x = [flipud(x_b(2:end, :)); x_f];
[status.closed, x] = geo.check_curve_loop(x, 'eps', h * 0.05);
status.fun_eval_cnt = status_f.fun_eval_cnt + status_b.fun_eval_cnt;
if ~status.closed
    status.completed = status_f.completed && status_b.completed;
else
    status.completed = true;
end
end

% ================================================================================
function [x, status] = search_direction(fdf, x0, direction, varargin)
% Helper function. Search a contour following one direction.

dim = size(x0, 2);

p = inputParser;
p.addParameter('method', 'euler', @(x) ischar(x));
p.addParameter('eps', 1e-8, @(x) validateattributes(x, {'double'}, {'scalar', 'positive'}));
p.addParameter('h', 0.05, @(x) validateattributes(x, {'double'}, {'scalar', 'positive'}));
p.addParameter('MaxPts', 100, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
p.parse(varargin{:});

fdx = @(x) wrap_to_fdx(fdf, x);
if strcmpi(p.Results.method, 'euler')
    method = @euler;
else
    error('Method name cannot be recognized!');
end

y0 = fdf(x0);

x = nan(p.Results.MaxPts, dim);
x(1, :) = x0;
h = p.Results.h;
h_min = h / 2^4;
h_max = h * 2^4;
idx = 1;
fun_eval_cnt = 1;

while idx < p.Results.MaxPts
    x0 = x(idx, :);

    x2_flag = false;
    while ~x2_flag && h > h_min
        [x1, x1_status] = method(fdx, x0, direction, h);
        fun_eval_cnt = fun_eval_cnt + x1_status.fun_eval_cnt;
        if isnan(x1_status.error)
            % Estimate the value error. Either from ODE method itself (like adaptive methods)
            % or compute it directly.
            y1 = fdf(x1);
            fun_eval_cnt = fun_eval_cnt + 1;
            x1_status.error = norm(y1 - y0);
        end

        [x2, x2_status] = ode.find_solution_fdf(fdf, x1, y0, 'eps', p.Results.eps);
        x2_flag = x2_status.finish;
        fun_eval_cnt = fun_eval_cnt + x2_status.fun_eval_cnt;

        % Then apply a simple adaptive scheme
        dx1 = norm(x2 - x1);
        if ~x2_flag || (dx1 > h * 0.05 || x1_status.error > max(p.Results.eps * 1e3, h * 0.3))
            h = h * 0.5;
        elseif dx1 < h * 0.01 && h < h_max * 0.5
            h = h * 2;
        end
    end

    if ~x2_flag
        break;
    else
        x(idx + 1, :) = x2;
        idx = idx + 1;
    end
end

[status.closed, x] = geo.check_curve_loop(x, 'eps', p.Results.h * 0.05);
status.completed = status.closed;
status.fun_eval_cnt = fun_eval_cnt;
if idx <= 1
    x = [];
elseif idx < p.Results.MaxPts
    status.completed = true;
end
end

% ================================================================================
function [y, dx] = wrap_to_fdx(fdf, x)
% Wrapper function. Convert fdf function to fdx function.
% An fdf function outputs value (y) and jacobian (jac), and an fdx function
% outputs value (y) and a contour tangent direction (dx)

[y, jac] = fdf(x);
dx = null(jac);

if det([jac', dx]) < 0
    dx = -dx;
end
dx = dx(:, 1)';
end

% ================================================================================
function [x1, s1] = euler(fdx, x0, direction, h)
% Simple Euler method.
% This method do NOT estimate error.

[~, dx] = fdx(x0);
x1 = x0 + direction * h * dx;
s1.error = nan;
s1.fun_eval_cnt = 1;
end
