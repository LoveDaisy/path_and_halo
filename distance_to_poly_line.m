function [dist, vn, idx_t] = distance_to_poly_line(qxy, xy)
% INPUT
%   qxy:    n*2
%   xy:     m*2

n = size(qxy, 1);
dim = size(qxy, 2);

dist = nan(n, 1);
vn = nan(n, dim);
idx_t = nan(n, 2);
for i = 1:n
    [tmp_dist, tmp_vn, tmp_t] = p2line_dist(qxy(i, :), xy(1:end-1, :), xy(2:end, :));
    [dist(i), tmp_idx] = min(tmp_dist);
    vn(i, :) = tmp_vn(tmp_idx, :);
    idx_t(i, :) = [tmp_idx, tmp_t(tmp_idx)];
end
end


function [dist, vp_n, t] = p2line_dist(p, xy1, xy2)
% INPUT
%   p:          1*2
%   xy1, xy2:   n*2

vp = bsxfun(@minus, p, xy1);
v2 = xy2 - xy1;

t = min(max(sum(vp .* v2, 2) ./ sum(v2.^2, 2), 0), 1);
vp_proj = bsxfun(@times, t, v2);
vp_n = bsxfun(@minus, vp, vp_proj);
dist = sqrt(sum(vp_n.^2, 2));
end