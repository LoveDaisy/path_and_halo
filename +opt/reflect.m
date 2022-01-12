function [ray_out, g] = reflect(ray_in, face_normal)
% INPUT
%  ray_in:       normalized, n*3
%  face_normal:  normalized
%  n0, n1:       refractive index

ray_n = size(ray_in, 1);

ray_out = nan(ray_n, 3);
g = nan(3, 3, ray_n);
valid_ind = sum(abs(ray_in), 2) > 1e-4;

if all(~valid_ind)
    return;
end

[ray_in, g_norm] = geo.normalize_vector(ray_in);

r = ray_in(valid_ind, :) - 2 * (ray_in(valid_ind, :) * face_normal') * face_normal;
g_r = (eye(3) - 2 * (face_normal' * face_normal)) * g_norm;

[r, g_norm_r] = geo.normalize_vector(r);
for i = 1:sum(valid_ind)
    g_r(:, :, i) = g_norm_r(:, :, i) * g_r(:, :, i);
end

ray_out(valid_ind, :) = r;
g(:, :, valid_ind) = g_r;
end
