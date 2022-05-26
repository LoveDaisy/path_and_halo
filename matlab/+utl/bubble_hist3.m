function bubble_hist3(data, grid_size, varargin)
% Plot a 3D bubble histogram
%
% INPUT
%   data:       m*3 or m*4, 3d data (with weights)
%   grid_size:  scalar or 1*3, 3d grid size
%
% OPTION
%   scale:      scalar, default 1.

p = inputParser;
p.addRequired('data', @(x) size(x, 2) == 3 || size(x, 2) == 4);
p.addRequired('grid_size', @(x) validateattributes(x, {'numeric'}, {'row'}));
p.addParameter('scale', 1, @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p.parse(data, grid_size, varargin{:});

[num, dim] = size(data);
if dim == 3
    w = ones(num, 1);
else
    w = data(:, 4);
end

bubble_xyz = bsxfun(@times, round(bsxfun(@times, data(:, 1:3), 1. / grid_size)), grid_size);
[bubble_xyz, ~, ic] = unique(bubble_xyz, 'rows');
bubble_s = accumarray(ic, w);

scatter3(bubble_xyz(:, 1), bubble_xyz(:, 3), bubble_xyz(:, 2), bubble_s * p.Results.scale + 1e-4);
end