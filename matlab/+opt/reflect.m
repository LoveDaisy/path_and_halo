function [ray_out, g] = reflect(ray_in, face_normal)
% INPUT
%  ray_in:       normalized, n*3
%  face_normal:  normalized
%  n0, n1:       refractive index

ray_n = size(ray_in, 1);
need_jacobian = nargout == 2;

ray_out = nan(ray_n, 3);
g = nan(3, 3, ray_n);
valid_ind = sum(abs(ray_in), 2) > 1e-4;
valid_cnt = sum(valid_ind);

if valid_cnt <= 0
    return;
end

r = ray_in(valid_ind, :) - 2 * (ray_in(valid_ind, :) * face_normal') * face_normal;
if need_jacobian
    g_r = nan(3, 3, valid_cnt);
    for i = 1:valid_cnt
        g_r(:, :, i) = (eye(3) - 2 * (face_normal' * face_normal));
    end
    g(:, :, valid_ind) = g_r;
end

ray_out(valid_ind, :) = r;
end
