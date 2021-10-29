function [ray_out, bending_angle, g_out, g_angle] = ...
    bending_angle_with_gradient(ray_in, face_normal, n)
% INPUT
%    ray_in:       m*2, [longitude, latitude] in degree
%    face_normal:  k*3, [nx, ny, nz] face normals
%    n:            k-vector, refract index after each face

n = [1; n(:)];
face_n = size(face_normal, 1);
ray_n = size(ray_in, 1);

% auto differentiation forward mode

% convert to xyz coordinate
ray_xyz = ll2xyz(ray_in);
g_xyz = ll2xyz_gradient(ray_in);

% refract on each face
curr_ray = ray_xyz;
g_ray = g_xyz;
for i = 1:face_n
    [curr_ray, g_curr] = refract_with_gradient(curr_ray, face_normal(i, :), n(i), n(i+1));
    for j = 1:ray_n
        g_ray(:, :, j) = g_curr(:, :, j) * g_ray(:, :, j);
    end
end

% convert back to longitude & latitude coordinate
ray_out = xyz2ll(curr_ray);
g_ll = xyz2ll_gradient(curr_ray);
g_out = zeros(2, 2, ray_n);
for i = 1:ray_n
    g_out(:, :, i) = g_ll(:, :, i) * g_ray(:, :, i);
end

% find out bending angle
cos_angle = sum(curr_ray .* ray_xyz, 2);
g_cos = zeros(ray_n, 2);
for i = 1:ray_n
    g_cos(i, :) = ray_xyz(i, :) * g_ray(:, :, i) + curr_ray(i, :) * g_xyz(:, :, i);
end
bending_angle = acosd(cos_angle);
g_angle = bsxfun(@times, -1 ./ sqrt(1 - cos_angle.^2), g_cos) * 180 / pi;
end


function xyz = ll2xyz(lonlat)
xyz = [cosd(lonlat(:, 2)) .* cosd(lonlat(:, 1)), ...
    cosd(lonlat(:, 2)) .* sind(lonlat(:, 1)), ...
    sind(lonlat(:, 2))];
end

function g = ll2xyz_gradient(lonlat)
% OUTPUT
%   g:      3*2*m
m = size(lonlat, 1);
g1 = [-cosd(lonlat(:, 2)) .* sind(lonlat(:, 1)), ...
    cosd(lonlat(:, 2)) .* cosd(lonlat(:, 1)), ...
    zeros(m, 1)];  % size m*3
g2 = [-sind(lonlat(:, 2)) .* cosd(lonlat(:, 1)), ...
    -sind(lonlat(:, 2)) .* sind(lonlat(:, 1)), ...
    cosd(lonlat(:, 2))];  % size m*3
g = cat(3, g1, g2) * pi / 180;  % size m*3*2
g = permute(g, [2, 3, 1]);
end


function ll = xyz2ll(xyz)
ll = [atan2d(xyz(:, 2), xyz(:, 1)), ...
    asind(xyz(:, 3) ./ sqrt(sum(xyz.^2, 2)))];
end

function g = xyz2ll_gradient(xyz)
% OUTPUT
%   g:      2*3*m
m = size(xyz, 1);
g1 = [-xyz(:, 2) ./ sum(xyz(:, 1:2).^2, 2), ...
    -xyz(:, 1) .* xyz(:, 3) ./ sqrt(sum(xyz(:, 1:2).^2, 2)) ./ sum(xyz.^2, 2)];  % size m*2
g2 = [xyz(:, 1) ./ sum(xyz(:, 1:2).^2, 2), ...
    -xyz(:, 2) .* xyz(:, 3) ./ sqrt(sum(xyz(:, 1:2).^2, 2)) ./ sum(xyz.^2, 2)];  % size m*2
g3 = [zeros(m, 1), ...
    sqrt(sum(xyz(:, 1:2).^2, 2)) ./ sum(xyz.^2, 2)];
g = cat(3, g1, g2, g3) * 180 / pi;  % size m*2*3
g = permute(g, [2, 3, 1]);
end
