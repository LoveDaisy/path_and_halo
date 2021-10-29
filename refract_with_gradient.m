function [ray_out, g] = refract_with_gradient(ray_in, face_normal, n)
% INPUT
%    ray_in:       m*2, [longitude, latitude] in degree
%    face_normal:  k*3, [nx, ny, nz] face normals
%    n:            k-vector, refract index after each face

n = [1; n(:)];
face_n = size(face_normal, 1);
ray_n = size(ray_in, 1);

% auto differentiation forward mode
dx = ones(size(ray_in));
x = ray_in;

% convert to xyz coordinate
dy = ll2xyz_gradient(x);
g = zeros(ray_n, 3);
for k = 1:ray_n
    g(k, :) =  dx(k, :) * dy(:, :, k)';
end
dx = g;
x = ll2xyz(x);

% refract on each face
for k = 1:face_n
    x = refract(x, face_normal(k, :), n(k), n(k+1));
end

% convert back to longitude & latitude coordinate
dy = xyz2ll_gradient(x);
g = zeros(ray_n, 2);
for k = 1:ray_n
    g(k, :) = dx(k, :) * dy(:, :, k)';
end
ray_out = xyz2ll(x);
end


function xyz = ll2xyz(lonlat)
xyz = [cosd(lonlat(:, 2)) .* cosd(lonlat(:, 1)), ...
    cosd(lonlat(:, 2)) .* sind(lonlat(:, 1)), ...
    sind(lonlat(:, 2))];
end

function g = ll2xyz_gradient(lonlat)
% OUTPUT
%   g:      3*2*k
m = size(lonlat, 1);
g1 = [-cosd(lonlat(:, 2)) .* sind(lonlat(:, 1)), ...
    cosd(lonlat(:, 2)) .* cosd(lonlat(:, 1)), ...
    zeros(m, 1)];  % size k*3
g2 = [-sind(lonlat(:, 2)) .* cosd(lonlat(:, 1)), ...
    -sind(lonlat(:, 2)) .* sind(lonlat(:, 1)), ...
    cosd(lonlat(:, 2))];  % size k*3
g = cat(3, g1, g2) * pi / 180;  % size k*3*2
g = permute(g, [2, 3, 1]);
end


function ll = xyz2ll(xyz)
ll = [atan2d(xyz(:, 2), xyz(:, 1)), ...
    asind(xyz(:, 3) ./ sqrt(sum(xyz.^2, 2)))];
end

function g = xyz2ll_gradient(xyz)
% OUTPUT
%   g:      2*3*k
m = size(xyz, 1);
g1 = [-xyz(:, 2) ./ sum(xyz(:, 1:2).^2, 2), ...
    -xyz(:, 1) .* xyz(:, 3) ./ sqrt(sum(xyz(:, 1:2).^2, 2)) ./ sum(xyz.^2, 2)];  % size k*2
g2 = [xyz(:, 1) ./ sum(xyz(:, 1:2).^2, 2), ...
    -xyz(:, 2) .* xyz(:, 3) ./ sqrt(sum(xyz(:, 1:2).^2, 2)) ./ sum(xyz.^2, 2)];  % size k*2
g3 = [zeros(m, 1), ...
    sqrt(sum(xyz(:, 1:2).^2, 2)) ./ sum(xyz.^2, 2)];
g = cat(3, g1, g2, g3) * 180 / pi;  % size k*2*3
g = permute(g, [2, 3, 1]);
end
