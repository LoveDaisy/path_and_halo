function [pts, s] = interp_curve(pts0, ds)
% Interpolate a curve with gievn interval ds
%
% INPUT
%   pts0:       n*d, end points for a polyline
%   ds:         scalar, minimumn length interval on curve
%
% OUTPUT
%   pts:        m*d, interpolated curve points
%   s:          m*1, length parameter for every interpolated point

len = sqrt(sum((pts0(1:end - 1, :) - pts0(2:end, :)).^2, 2));
s0 = cumsum([0; len]);

len_diff = inf;
while len_diff > 0.1 * ds
    s = unique([(0:ds:s0(end))'; s0], 'sorted');
    pts = interp1(s0, pts0, s, 'spline');

    len = sqrt(sum((pts(1:end - 1, :) - pts(2:end, :)).^2, 2));
    tmp_s = cumsum([0; len]);
    tmp_d = pdist2(pts0, pts);
    [~, tmp_idx] = min(tmp_d, [], 2);
    len_diff = abs(tmp_s(end) - s(end));
    
    s0 = tmp_s(tmp_idx);
end
s = tmp_s;
end
