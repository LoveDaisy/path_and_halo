function [vec, g_vec] = normalize_vector_with_gradient(vec)
% INPUT
%   vec:        n*d
% OUTPUT
%   vec:        n*d
%   g_vec:      d*d*n

num = size(vec, 1);

sq_vec = vec.^2;
sq_vec_norm = sum(sq_vec, 2);
rcp_vec_norm = 1 ./ sqrt(sq_vec_norm);

g_vec = [sq_vec_norm - sq_vec(:, 1), -vec(:, 1) .* vec(:, 2), -vec(:, 1) .* vec(:, 3), ...
    -vec(:, 2) .* vec(:, 1), sq_vec_norm - sq_vec(:, 2), -vec(:, 2) .* vec(:, 3), ...
    -vec(:, 3) .* vec(:, 1), -vec(:, 3) .* vec(:, 2), sq_vec_norm - sq_vec(:, 3)];
g_vec = bsxfun(@times, g_vec, rcp_vec_norm.^3)';
g_vec = reshape(g_vec, [3, 3, num]);

vec = bsxfun(@times, vec, rcp_vec_norm);
end