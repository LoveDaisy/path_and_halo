function [ray_out, bending_angle, g_out, g_angle] = ...
    trace_ray_xyz_with_gradient(ray_in, face_normal, n)
% INPUT
%    ray_in:       m*3, [longitude, latitude] in degree
%    face_normal:  k*3, [nx, ny, nz] face normals
%    n:            k-vector, refract index after each face

n = [1; n(:)];
face_n = size(face_normal, 1);
ray_n = size(ray_in, 1);

% normalize
ray_in_norm = sqrt(sum(ray_in.^2, 2));
g_norm = [sum(ray_in(:, [2, 3]).^2, 2), -prod(ray_in(:, [1, 2]), 2), -prod(ray_in(:, [1, 3]), 2), ...
    -prod(ray_in(:, [2, 1]), 2), sum(ray_in(:, [1, 3]).^2, 2), -prod(ray_in(:, [2, 3]), 2), ...
    -prod(ray_in(:, [3, 1]), 2), -prod(ray_in(:, [3, 2]), 2), sum(ray_in(:, [1, 2]).^2, 2)];
g_norm = bsxfun(@times, g_norm, 1 ./ ray_in_norm.^3);
g_norm = reshape(g_norm', [3, 3, ray_n]);
ray_in = bsxfun(@times, ray_in, 1 ./ ray_in_norm);

% refract on each face
curr_ray = ray_in;
g_out = g_norm;
for i = 1:face_n
    [curr_ray, g_curr] = refract_with_gradient(curr_ray, face_normal(i, :), n(i), n(i+1));
    for j = 1:ray_n
        g_out(:, :, j) = g_curr(:, :, j) * g_out(:, :, j);
    end
end
ray_out = curr_ray;

% find out bending angle
cos_angle = sum(curr_ray .* ray_in, 2);
g_cos = zeros(ray_n, 3);
for i = 1:ray_n
    g_cos(i, :) = ray_in(i, :) * g_out(:, :, i) + curr_ray(i, :) * g_norm(:, :, i);
end
bending_angle = acosd(cos_angle);
g_angle = bsxfun(@times, -1 ./ sqrt(1 - cos_angle.^2), g_cos) * 180 / pi;
end

