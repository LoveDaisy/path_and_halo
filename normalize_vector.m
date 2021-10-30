function v = normalize_vector(v)
% INPUT
%   v:      n*m, m-demension vectors

v = bsxfun(@times, v, 1./sqrt(sum(v.^2, 2)));
end