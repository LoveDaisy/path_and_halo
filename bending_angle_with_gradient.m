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
[ray_xyz, g_xyz] = ll2xyz_with_gradient(ray_in);

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
[ray_out, g_ll] = xyz2ll_with_gradient(curr_ray);
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

