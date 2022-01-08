function [w, interp_s, interp_p, interp_rot] = compute_llr_weight(rot_contour, axis_pdf, config)
[p0, det_j0, face_factor0, t_factor0] = compute_components(rot_contour, axis_pdf, config);
interest_idx = face_factor0 >= 1e-10;
idx1 = find(diff(double(interest_idx)) > 0);
idx2 = find(diff(double(interest_idx)) < 0) + 1;
if interest_idx(1)
    idx1 = [1; idx1];
end
if interest_idx(end)
    idx2 = [idx2; length(p0)];
end

w = 0;
interp_s = [];
interp_p = [];
interp_rot = [];
if sum(~isnan(det_j0)) < 2
    return;
end

[interp_rot, s, interp_s] = spline_interp_rot(rot_contour);
if isempty(interp_s)
    return;
end
interp_p = axis_pdf(interp_rot);
interp_det_j = exp(interp1(s, log(det_j0), interp_s, 'linear', 'extrap'));
interp_face_factor = interp1(s, face_factor0, interp_s, 'linear', 'extrap');
interp_t_factor = exp(interp1(s, log(t_factor0), interp_s, 'linear', 'extrap'));
for i = 1:length(idx1)
    i1 = find(interp_s <= s(idx1(i)), 1, 'last');
    i2 = find(interp_s >= s(idx2(i)), 1, 'first');
    if isempty(i1)
        i1 = 1;
    end
    if isempty(i2)
        i2 = length(interp_s);
    end

    tmp_idx = false(size(interp_p));
    tmp_idx(i1:i2) = true;
    tmp_idx = tmp_idx & interp_p > 1e-10;
    if sum(tmp_idx) > 100 && (mean(face_factor0(idx1:idx2)) > 2e-1 || ...
            mean(interp_face_factor(i1:i2) .* interp_p(i1:i2)) < 1e-8 || ...
            abs(interp_s(i1) - interp_s(i2)) > 50)
        continue;
    end

    [~, tmp_det_j, tmp_face_factor, tmp_t_factor] = compute_components(interp_rot(tmp_idx, :), axis_pdf, config);
    interp_det_j(tmp_idx) = tmp_det_j;
    interp_face_factor(tmp_idx) = tmp_face_factor;
    interp_t_factor(tmp_idx) = tmp_t_factor;
end

interp_p = [interp_p ./ interp_det_j .* interp_face_factor .* interp_t_factor, ...
            interp_p, 1 ./ interp_det_j, interp_face_factor, interp_t_factor];
w = sum((interp_p(1:end - 1, 1) + interp_p(2:end, 1)) / 2 .* diff(interp_s));
end

function [p, det_j, face_factor, t_factor] = compute_components(rot, axis_pdf, config)
sun_ll = config.sun_ll;
ray_in_ll = [sun_ll(1) + 180, -sun_ll(2)];
ray_in_xyz = geo.ll2xyz(ray_in_ll);

p = axis_pdf(rot);

det_j = zeros(size(p));
face_factor = zeros(size(p));
t_factor = zeros(size(p));

for i = 1:length(p)
    [~, curr_jacob] = opt.crystal_system(rot(i, :), ray_in_ll, config.crystal, config.trace);

    det_j(i) = jacobian_factor(curr_jacob(:, :, 1));
    tmp_rot_mat = rotz(90 + rot(i, 1)) * rotx(90 - rot(i, 2)) * rotz(rot(i, 3));

    tmp_crystal = config.crystal;
    tmp_crystal.face_norm = config.crystal.face_norm * tmp_rot_mat';
    tmp_crystal.vtx = config.crystal.vtx * tmp_rot_mat';

    tmp_face_area = entry_face_factor(tmp_crystal, ray_in_xyz);
    face_factor(i) = tmp_face_area(config.trace.fid(1));

    [t, a] = transit_face_factor(ray_in_xyz, tmp_crystal, config.trace);
    t_factor(i) = t;
    face_factor(i) = face_factor(i) * a;
end
end

function j = jacobian_factor(j_rot)
m = j_rot(:, :) * j_rot(:, :)';
a0 = j_rot(1, :);
a = a0 / norm(a0);

b = a * j_rot(:, :)';
b = [b(2), -b(1)];
b = b / sqrt(b * m * b');

min_det = 1e-4;
j = abs(det(m) * det([[1 / norm(a0); 0], b'])) + min_det;
if isnan(j)
    j = min_det;
end
end

function a = entry_face_factor(crystal, ray_in_xyz)
face_area = crystal.face_area .* (-crystal.face_norm * ray_in_xyz');
a = max(face_area, 0) / sum(max(face_area, 0));
end

function [t_factor, face_factor] = transit_face_factor(ray_in_xyz, crystal, trace)
trace_n = opt.generate_trace_n(crystal, trace);
vtx0 = crystal.vtx(crystal.face{trace.fid(1)}, :);
r = ray_in_xyz;

t_factor = 1;
face_factor = 1;
for i = 1:length(trace.fid)
    n0 = trace_n(i);
    n1 = trace_n(i + 1);
    fn = crystal.face_norm(trace.fid(i), :);
    if n0 * n1 > 0
        [t, r, q] = transit_factor(r, fn, n0, n1);
    elseif abs(n0) >= 1 - 1e-6
        [t, ~, ~] = transit_factor(r, fn, abs(n0), 1);
        if isnan(t)
            t = 1;
        else
            t = 1 - t;
        end
        r = opt.reflect(r, fn);
        q = 1;
    else
        [t, ~, ~] = transit_factor(r, fn, 1, crystal.n);
        t = 1 - t;
        q = 1;
        r = opt.reflect(r, fn);
    end
    t_factor = t_factor * t * q;

    if i < length(trace.fid) && face_factor > 1e-8
        vtx1 = crystal.vtx(crystal.face{trace.fid(i + 1)}, :);
        [a, tmp_vtx] = face_intersection_factor(vtx0, vtx1, r);
        face_factor = face_factor * a;
        vtx0 = tmp_vtx;
    end
end
end

function [t, ray_out, q_factor] = transit_factor(ray_in, face_normal, n0, n1)
if n0 * n1 > 0
    ray_out = opt.refract(ray_in, face_normal, n0, n1);
else
    ray_out = opt.reflect(ray_in, face_normal);
end
n0 = abs(n0);
n1 = abs(n1);
cos_qi = abs(dot(ray_in, face_normal));
cos_qt = abs(dot(ray_out, face_normal));
Rs1 = abs((n0 * cos_qi - n1 * cos_qt) / (n0 * cos_qi + n1 * cos_qt))^2;
Rp1 = abs((n0 * cos_qt - n1 * cos_qi) / (n0 * cos_qt + n1 * cos_qi))^2;
t = (1 - (Rs1 + Rp1) / 2);
q_factor = cos_qi / cos_qt;
end

function [area_factor, vtx] = face_intersection_factor(vtx_q, vtx_0, ray_xyz)
% INPUT
%   vtx_q:      entry face vertex
%   vtx_0:      exit face vertex. it must be *CONVEX*
%   ray_xyz:    projection ray

e1 = vtx_0(2, :) - vtx_0(1, :);
e2 = vtx_0(4, :) - vtx_0(1, :);
e3 = cross(e1, e2);

uv0 = bsxfun(@minus, vtx_0, vtx_0(1, :)) / [e1; e2; e3];
uv0 = uv0(:, 1:2);
uvt = bsxfun(@minus, vtx_q, vtx_0(1, :)) / [e1; e2; -ray_xyz];
uvq = uvt(:, 1:2);

mask_vtx_num = size(uv0, 1);
q_vtx_num = size(uvq, 1);
uv = nan(q_vtx_num * 2, 2);

z_eps = 1e-8;
for i1 = 1:mask_vtx_num
    i2 = mod(i1, mask_vtx_num) + 1;
    i3 = mod(i2, mask_vtx_num) + 1;

    z0 = (uv0(i2, 1) - uv0(i1, 1)) * (uv0(i3, 2) - uv0(i2, 2)) - ...
        (uv0(i2, 2) - uv0(i1, 2)) * (uv0(i3, 1) - uv0(i2, 1));
    z2 = nan;
    uv_i = 1;
    for j1 = 1:q_vtx_num
        j2 = mod(j1, q_vtx_num) + 1;
        if ~isnan(z2)
            z1 = z2;
        else
            z1 = (uv0(i2, 1) - uv0(i1, 1)) * (uvq(j1, 2) - uv0(i1, 2)) - ...
                (uv0(i2, 2) - uv0(i1, 2)) * (uvq(j1, 1) - uv0(i1, 1));
        end
        z2 = (uv0(i2, 1) - uv0(i1, 1)) * (uvq(j2, 2) - uv0(i1, 2)) - ...
            (uv0(i2, 2) - uv0(i1, 2)) * (uvq(j2, 1) - uv0(i1, 1));

        if z0 * z1 > z_eps && z0 * z2 < z_eps
            % From inner to outer: record tow points
            uv(uv_i, :) = uvq(j1, :);
            ab =- (uv0(i1, :) - uvq(j1, :)) / [uv0(i2, :) - uv0(i1, :); uvq(j2, :) - uvq(j1, :)];
            uv(uv_i + 1, :) = uv0(i1, :) + ab(1) * (uv0(i2, :) - uv0(i1, :));
            uv_i = uv_i + 2;
        elseif z0 * z1 < -z_eps && z0 * z2 > -z_eps
            % From outer to inner: record one point
            ab =- (uv0(i1, :) - uvq(j1, :)) / [uv0(i2, :) - uv0(i1, :); uvq(j2, :) - uvq(j1, :)];
            uv(uv_i, :) = uv0(i1, :) + ab(1) * (uv0(i2, :) - uv0(i1, :));
            uv_i = uv_i + 1;
        elseif z0 * z1 > z_eps && z0 * z2 > z_eps
            % Two inner: record one point
            uv(uv_i, :) = uvq(j1, :);
            uv_i = uv_i + 1;
        elseif abs(z1) < z_eps
            uv(uv_i, :) = uvq(j1, :);
            uv_i = uv_i + 1;
        end
    end
    q_vtx_num = uv_i - 1;
    uvq = uv(1:q_vtx_num, :);
end
uv = uv(1:q_vtx_num, :);
vtx = uv * [vtx_0(2, :) - vtx_0(1, :); vtx_0(4, :) - vtx_0(1, :)];
vtx = bsxfun(@plus, vtx, vtx_0(1, :));

if size(uv, 1) > 2
    area1 = sum(uv(:, 1) .* [uv(end, 2); uv(1:end - 1, 2)]) - ...
        sum([uv(end, 1); uv(1:end - 1, 1)] .* uv(:, 2));
    area0 = sum(uvt(:, 1) .* [uvt(end, 2); uvt(1:end - 1, 2)]) - ...
        sum([uvt(end, 1); uvt(1:end - 1, 1)] .* uvt(:, 2));
    area_factor = area1 / area0;
else
    area_factor = 0;
end
end

function [interp_rot, s, interp_s] = spline_interp_rot(rot_contour)
ds = sqrt(sum(diff(rot_contour).^2, 2));
s = [0; cumsum(ds)];
num = size(rot_contour, 1);
ss = 0.05;

if max(s) < ss
    interp_rot = [];
    interp_s = [];
    return;
end

is_period = norm(rot_contour(1, :) - rot_contour(end, :)) < 1e-8;

s_diff = inf;
iter_num = 1;
while s_diff > 1e-3 && iter_num < 3
    interp_s = [(s(1):ss:s(end))'; s];
    [c, ~, ic] = unique(interp_s);
    [interp_s, ics] = sort(c);
    [~, isc] = sort(ics);
    if is_period
        interp_rot = interp1([s; s(end) + s(2:end)], [rot_contour; rot_contour(2:end, :)], interp_s, 'spline');
        interp_rot = interp_rot(1:length(interp_s), :);
    else
        interp_rot = interp1(s, rot_contour, interp_s, 'spline');
    end

    new_ds = sqrt(sum(diff(interp_rot).^2, 2));
    new_s = [0; cumsum(new_ds)];
    s_diff = abs(new_s(end) - s(end)) / s(end);

    s = new_s(isc(ic(end - num + 1:end)));
    iter_num = iter_num + 1;
end
end
