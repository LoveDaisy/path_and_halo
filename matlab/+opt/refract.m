function [ray_out, g] = refract(ray_in, face_normal, n0, n1)
% INPUT
%  ray_in:       normalized, n*3
%  face_normal:  normalized
%  n0, n1:       refractive index

ray_n = size(ray_in, 1);
need_jacobian = nargout == 2;

ray_out = nan(ray_n, 3);
g = nan(3, 3, ray_n);
valid_ind = sum(abs(ray_in), 2) > 1e-4;

% [ray_in, g_norm] = geo.normalize_vector(ray_in);

cos_alpha = ray_in * face_normal';
g_cos_alpha = nan(ray_n, 3);
if need_jacobian
    for i = 1:ray_n
        g_cos_alpha(i, :) = face_normal';
    end
end

delta = cos_alpha.^2 - 1 + (n1 / n0)^2;
g_delta = g_cos_alpha;
if need_jacobian
    for i = 1:ray_n
        g_delta(i, :) = g_delta(i, :) * cos_alpha(i) * 2;
    end
end

valid_delta = delta > 0;
valid_ind = valid_ind & valid_delta;
g_delta(~valid_ind, :) = nan;
g_cos_alpha(~valid_ind, :) = nan;

if ~any(valid_ind)
    return;
end

valid_ind = find(valid_ind);
a = cos_alpha(valid_ind) - sign(cos_alpha(valid_ind)) .* sqrt(delta(valid_ind));
g_a = g_cos_alpha(valid_ind, :);
if need_jacobian
    for i = 1:length(valid_ind)
        g_a(i, :) = g_a(i, :) - sign(cos_alpha(valid_ind(i))) / 2 / sqrt(delta(valid_ind(i))) * ...
            g_delta(valid_ind(i), :);
    end
end

r = ray_in(valid_ind, :) - a * face_normal;
g_r = nan(3, 3, length(valid_ind));
if need_jacobian
    for i = 1:length(valid_ind)
        g_r(:, :, i) = eye(3) - face_normal' * g_a(i, :);
    end

    [r, g_norm_r] = geo.normalize_vector(r);
    for i = 1:length(valid_ind)
        g_r(:, :, i) = g_norm_r(:, :, i) * g_r(:, :, i);
    end
else
    r = geo.normalize_vector(r);
end

ray_out(valid_ind, :) = r;
g(:, :, valid_ind) = g_r;
end
