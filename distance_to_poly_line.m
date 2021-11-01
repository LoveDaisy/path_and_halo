function d = distance_to_poly_line(qxy, xy)
% INPUT
%   qxy:    n*2
%   xy:     m*2

n = size(qxy, 1);
d = nan(n, 1);
for i = 1:n
    d(i) = min(p2line_dist(qxy(i, :), xy(1:end-1, :), xy(2:end, :)));
end
end


function d = p2line_dist(p, xy1, xy2)
% INPUT
%   p:          1*2
%   xy1, xy2:   n*2

vp = bsxfun(@minus, p, xy1);
v2 = xy2 - xy1;

vp_proj = min(max(sum(vp .* v2, 2) ./ sum(v2.^2, 2), 0), 1);
vp_proj = bsxfun(@times, vp_proj, v2);
d = sqrt(sum(bsxfun(@minus, vp, vp_proj).^2, 2));
end