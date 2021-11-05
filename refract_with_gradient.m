function [ray_out, g] = refract_with_gradient(ray_in, face_normal, n0, n1)
% INPUT
%  ray_in:       normalized, n*3
%  face_normal:  normalized
%  n0, n1:       refractive index

ray_n = size(ray_in, 1);

ray_out = nan(size(ray_in));
g = nan(3, 3, ray_n);
ray_in_norm = sqrt(sum(ray_in.^2, 2));
valid_input = sum(abs(ray_in), 2) > 1e-4;

cos_alpha = ray_in * face_normal' ./ ray_in_norm;
g_cos_alpha = [[sum(ray_in(:, [2, 3]).^2, 2), -prod(ray_in(:, [1, 2]), 2), ...
    -prod(ray_in(:, [1, 3]), 2)] * face_normal', ...
    [-prod(ray_in(:, [2, 1]), 2), sum(ray_in(:, [1, 3]).^2, 2), -prod(ray_in(:, [2, 3]), 2)] * ...
    face_normal', ...
    [-prod(ray_in(:, [3, 1]), 2), -prod(ray_in(:, [3, 2]), 2), sum(ray_in(:, [1, 2]).^2, 2)] * ...
    face_normal'];
g_cos_alpha = bsxfun(@times, g_cos_alpha, 1 ./ ray_in_norm.^3);

delta = cos_alpha .^ 2 - 1 + (n1 / n0) ^ 2;
g_delta = bsxfun(@times, cos_alpha * 2, g_cos_alpha);

valid_delta = delta > 0;
valid_ind = valid_input & valid_delta;
g_delta(~valid_ind, :) = nan;
g_cos_alpha(~valid_ind, :) = nan;

if ~any(valid_ind)
    return;
end
a = cos_alpha(valid_ind) - sign(cos_alpha(valid_ind)) .* sqrt(delta(valid_ind));
g_a = g_cos_alpha(valid_ind, :) - bsxfun(@times, sign(cos_alpha(valid_ind)) / ...
    2 ./ sqrt(delta(valid_ind)), g_delta(valid_ind, :));

r = bsxfun(@times, ray_in(valid_ind, :), 1 ./ ray_in_norm(valid_ind)) - bsxfun(@times, a, face_normal);
g_r1 = bsxfun(@times, [sum(ray_in(valid_ind, [2, 3]).^2, 2), ...
    -prod(ray_in(valid_ind, [1, 2]), 2), ...
    -prod(ray_in(valid_ind, [1, 3]), 2)], 1 ./ ray_in_norm(valid_ind).^3) - face_normal(1) * g_a;
g_r2 = bsxfun(@times, [-prod(ray_in(valid_ind, [1, 2]), 2), ...
    sum(ray_in(valid_ind, [1, 3]).^2, 2), ...
    -prod(ray_in(valid_ind, [2, 3]), 2)], 1 ./ ray_in_norm(valid_ind).^3) - face_normal(2) * g_a;
g_r3 = bsxfun(@times, [-prod(ray_in(valid_ind, [1, 3]), 2), ...
    -prod(ray_in(valid_ind, [2, 3]), 2), ...
    sum(ray_in(valid_ind, [1, 2]).^2, 2)], 1 ./ ray_in_norm(valid_ind).^3) - face_normal(3) * g_a;
g_r = cat(3, g_r1, g_r2, g_r3);
g_r0 = permute(g_r, [3, 2, 1]);

r_norm = sqrt(sum(r.^2, 2));
g_r1 = bsxfun(@times, [sum(r(:, [2, 3]).^2, 2), ...
    -prod(r(:, [1, 2]), 2), ...
    -prod(r(:, [1, 3]), 2)], 1 ./ r_norm.^3);
g_r2 = bsxfun(@times, [-prod(r(:, [2, 1]), 2), ...
    sum(r(:, [1, 3]).^2, 2), ...
    -prod(r(:, [2, 3]), 2)], 1 ./ r_norm.^3);
g_r3 = bsxfun(@times, [-prod(r(:, [3, 1]), 2), ...
    -prod(r(:, [3, 2]), 2), ...
    sum(r(:, [1, 2]).^2, 2)], 1 ./ r_norm.^3);  % n*3
g_r = cat(3, g_r1, g_r2, g_r3);
g_r = permute(g_r, [3, 2, 1]);
r = bsxfun(@times, r, 1 ./ r_norm);
for i = 1:sum(valid_ind)
    g_r(:, :, i) = g_r(:, :, i) * g_r0(:, :, i);
end

ray_out(valid_ind, :) = r;
g(:, :, valid_ind) = g_r;
end