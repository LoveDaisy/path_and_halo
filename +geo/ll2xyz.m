function [xyz, jac] = ll2xyz(lon_lat)
% Convert spherical coordinate [longitude, latitude] into Cartesian coordinate [x, y, z]
% Clearly [x, y, z] lies on the unit sphere.
%
% INPUT
%   lon_lat:        n*2, [longitude, latitude]
%
% OUTPUT
%   xyz:            n*3, [x, y, z]
%   jac:            3*2*n, Jacobian of this transformation

p = inputParser;
p.addRequired('xyz', @(x) validateattributes(x, {'numeric'}, {'2d', 'ncols', 2}));
p.parse(lon_lat);

xyz = [cosd(lon_lat(:, 2)) .* cosd(lon_lat(:, 1)), ...
        cosd(lon_lat(:, 2)) .* sind(lon_lat(:, 1)), ...
        sind(lon_lat(:, 2))];

m = size(lon_lat, 1);
g1 = [-cosd(lon_lat(:, 2)) .* sind(lon_lat(:, 1)), ...
        cosd(lon_lat(:, 2)) .* cosd(lon_lat(:, 1)), ...
        zeros(m, 1)]; % size m*3
g2 = [-sind(lon_lat(:, 2)) .* cosd(lon_lat(:, 1)), ...
        -sind(lon_lat(:, 2)) .* sind(lon_lat(:, 1)), ...
        cosd(lon_lat(:, 2))]; % size m*3
jac = cat(3, g1, g2) * pi / 180; % size m*3*2
jac = permute(jac, [2, 3, 1]);
end
