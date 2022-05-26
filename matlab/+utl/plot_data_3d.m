function plot_data_3d(data, sub_idx, varargin)
% A simple wrapper to plot multi-dimention data in a 3d subspace.

if isempty(sub_idx)
    sub_idx = [1, 2, 3];
end
plot3(data(:, sub_idx(1)), data(:, sub_idx(2)), data(:, sub_idx(3)), varargin{:});
end