function [vec, g_vec] = normalize_vector(vec)
% INPUT
%   vec:        n*d
% OUTPUT
%   vec:        n*d
%   g_vec:      d*d*n

num = size(vec, 1);
d = size(vec, 2);
need_jacobian = nargout == 2;

sq_vec = vec.^2;
sq_vec_norm = sum(sq_vec, 2);
rcp_vec_norm = 1 ./ sqrt(sq_vec_norm);

g_vec = nan(d, d, num);
for i = 1:num
    if need_jacobian
        tmp_g = -vec(i, :)' * vec(i, :);
        tmp_g = tmp_g * rcp_vec_norm(i)^3;
        for j = 1:d
            tmp_g(j, j) = tmp_g(j, j) + rcp_vec_norm(i);
        end
        g_vec(:, :, i) = tmp_g;
    end

    vec(i, :) = vec(i, :) * rcp_vec_norm(i);
end
end
