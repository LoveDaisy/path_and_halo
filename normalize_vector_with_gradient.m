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

for i = 1:num
    g_vec(i, :) = g_vec(i, :) * rcp_vec_norm(i)^3;
    vec(i, :) = vec(i, :) * rcp_vec_norm(i);
end
g_vec = reshape(g_vec', [3, 3, num]);
end