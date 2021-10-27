function ray_out = refract(ray_in, face_normal, n0, n1)
% INPUT
%  ray_in:       normalized, n*3
%  face_normal:  normalized
%  n0, n1:       refractive index

ray_out = zeros(size(ray_in));
valid_input = sum(abs(ray_in), 2) > 1e-4;
if ~any(valid_input)
    return;
end

cos_alpha = ray_in * face_normal';
delta = cos_alpha .^ 2 - 1 + (n1 / n0) ^ 2;
valid_delta = delta > 0;
valid_ind = valid_input & valid_delta;

if ~any(valid_ind)
    return;
end
a = cos_alpha(valid_ind) - sign(cos_alpha(valid_ind)) .* sqrt(delta(valid_ind));
r = ray_in(valid_ind, :) - bsxfun(@times, a, face_normal);
r = bsxfun(@times, r, 1 ./ sqrt(sum(r.^2, 2)));
ray_out(valid_ind, :) = r;
end