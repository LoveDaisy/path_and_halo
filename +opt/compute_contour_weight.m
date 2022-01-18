function [weight, cmp_interp, rot_interp] = compute_contour_weight(rot_contour, axis_pdf, config)
% Compute the total weight of a rotation contour.
% Weight components are: [axis_prob, det_jac, geo_factor, transit_factor]. See help of
% compute_weight_components() below for detail.
%
% INPUT
%   rot_contour:        n*3 or n*4
%   axis_pdf:           function handle, take LLR as input and return probability
%   config:             struct, with fields of:
%       .crystal        struct
%       .trace          struct
%       .sun_ll         sun position, [longitude, latitude] in degree
%
% OUTPUT
%   weight:             scalar
%   cmp:                m*6, [s, w, axis_p, det_jac, geo_factor, transit_factor]

dim = size(rot_contour, 2);

cmp0 = compute_weight_components(rot_contour, axis_pdf, config);

% Calculate length of polyline as initial value
len = sum(sqrt(sum(diff(rot_contour).^2, 2)));

% Interpolate rotations
[rot_interp, s_interp, s0, interp_idx] = geo.interp_curve(rot_contour, len * 0.01);
interp_num = size(rot_interp, 1);

% Convert to LLR space
if dim == 4
    rot_interp = geo.quat2llr(rot_interp);
    diff_s = sqrt(sum(diff(rot_interp).^2, 2));
    s_interp = [0; cumsum(diff_s)];
    discontinuity_idx = find(diff_s > 30);
    discontinuity_idx = [0; discontinuity_idx; interp_num];
    for i = 2:length(discontinuity_idx)-1
        idx1 = discontinuity_idx(i)+1;
        idx2 = discontinuity_idx(i+1);
        s_interp(idx1:idx2) = s_interp(idx1:idx2) - diff_s(idx1 - 1) + diff_s(idx1 - 2);
    end
    s0 = s_interp(interp_idx);
elseif dim ~= 3
    error('Input rotation must have dimesion of 3 or 4!');
end

if norm(rot_contour(1, :) - rot_contour(end, :)) < 1e-10
    % For closed loop
    s0(end) = s_interp(end);
end

% Interpolate components as initial value
cmp_interp = nan(length(s_interp), 6);
cmp_interp(:, 1) = s_interp;
cmp_interp(:, 3) = axis_pdf(rot_interp);
cmp_interp(:, 4) = exp(interp1(s0, log(cmp0(:, 2)), s_interp, 'linear', 'extrap'));
cmp_interp(:, 5) = interp1(s0, cmp0(:, 3), s_interp, 'linear', 'extrap');
cmp_interp(:, 6) = exp(interp1(s0, log(cmp0(:, 4)), s_interp, 'linear', 'extrap'));

% Find out those with large probability && small geometry factor && small transit_factor
interest_idx = cmp_interp(:, 3) >= 1e-12 & (cmp_interp(:, 5) <= 1e-1 | cmp_interp(:, 6) <= 1e-1);
idx1 = find(diff(double(interest_idx)) > 0);
idx2 = find(diff(double(interest_idx)) < 0) + 1;
if interest_idx(1)
    idx1 = [1; idx1];
end
if interest_idx(end)
    idx2 = [idx2; interp_num];
end

% Compute components only on interested segments
for i = 1:length(idx1)
    i1 = idx1(i);
    i2 = idx2(i);
    cmp_interp(i1:i2, 3:6) = compute_weight_components(rot_interp(i1:i2, :), axis_pdf, config);
end
cmp_interp = max(cmp_interp, 0);
cmp_interp(:, 2) = cmp_interp(:, 3) .* cmp_interp(:, 4) .* cmp_interp(:, 5) .* cmp_interp(:, 6);
weight = sum(diff(s_interp) .* (cmp_interp(1:end - 1, 2) + cmp_interp(2:end, 2)) / 2);
end

function cmp = compute_weight_components(rot, axis_pdf, config)
% Compute all weight components. These factors are:
% 1. Axis distribution probability
% 2. Determinant of Jacobian.
% 3. Geometry factor (face area factor on entry and during ray propagation)
% 4. Ray transit factor (refract and reflect loss, following Fresnel's equation)
%
% INPUT
%   rot:            n*3 or n*4
%   axis_pdf:       function handle, take LLR as input and return probability
%   config:         struct, with fields of:
%       .crystal    struct
%       .trace      struct
%       .sun_ll     sun position, [longitude, latitude] in degree
%
% OUTPUT
%   cmp:            n*4, [axis_p, det_jac, geo_factor, transit_factor]

rot_dim = size(rot, 2);
if rot_dim == 3
    llr = geo.normalize_llr(rot);
elseif rot_dim == 4
    llr = geo.quat2llr(rot);
else
    error('Rotations must have dimenstion 3 or 4!');
end

axis_prob = axis_pdf(llr);
det_jac = compute_jacobian_deternimant(llr, config);
entry_factor = entry_face_factor(llr, config);
[transit_factor, geo_factor] = transit_geo_factor(llr, config);

cmp = [axis_prob, 1 ./ det_jac, entry_factor .* geo_factor, transit_factor];
end

% ================================================================================
function det_jac = compute_jacobian_deternimant(llr, config)
% Compute deteminant of Jacobian. For LLR representation, Jacobian is 2*2 matrix.
% What does this Jacobian mean? We simply image this, when a point moves within a small
% region perpendicular to contour (two degree of freedom), the output will also moves with in
% a small region (two degree of freedom), and thus Jacobian is 2*2 matrix.
% And clearly, the determinant of Jacobian, is the ratio of these two reagion area.

rot_num = size(llr, 1);
ray_in_ll = [config.sun_ll(1) + 180, -config.sun_ll(2)];
fdf = @(rot) opt.crystal_system(rot, ray_in_ll, config.crystal, config.trace);

[~, jac] = fdf(llr);
det_jac = nan(size(llr, 1), 1);
for i = 1:rot_num
    curr_jac = jac(:, :, i);

    m = curr_jac(:, :) * curr_jac(:, :)';
    a0 = curr_jac(1, :);
    a = a0 / norm(a0);

    b = a * curr_jac(:, :)';
    b = [b(2), -b(1)];
    b = b / sqrt(b * m * b');

    min_det = 1e-4;
    j = abs(det(m) * det([[1 / norm(a0); 0], b'])) + min_det;
    if isnan(j)
        j = min_det;
    end
    det_jac(i) = j;
end
end

% ================================================================================
function entry_factor = entry_face_factor(llr, config)
num = size(llr, 1);
rot_mat = geo.llr2mat(llr);

sun_ll = config.sun_ll;
ray_in_ll = [sun_ll(1) + 180, -sun_ll(2)];
ray_in_xyz = geo.ll2xyz(ray_in_ll);

entry_factor = nan(num, 1);
entry_fid = config.trace.fid(1);
for i = 1:num
    face_norm = config.crystal.face_norm * rot_mat(:, :, i)';
    face_area = config.crystal.face_area .* (-face_norm * ray_in_xyz');
    entry_factor(i) = max(face_area(entry_fid), 0) / sum(max(face_area, 0));
end
end

% ================================================================================
function [transit_factor, geo_factor] = transit_geo_factor(llr, config)
crystal = config.crystal;
trace = config.trace;
trace_n = opt.generate_trace_n(crystal, trace);

ray_in_ll = [config.sun_ll(1) + 180, -config.sun_ll(2)];
ray_in_xyz = geo.ll2xyz(ray_in_ll);

face_cnt = length(trace_n) - 1;
num = size(llr, 1);
rot_mat = geo.llr2mat(llr);

transit_factor = ones(num, 1);
geo_factor = ones(num, 1);

for i = 1:num
    r = ray_in_xyz;
    face_norm = crystal.face_norm * rot_mat(:, :, i)';
    crystal_vtx = crystal.vtx * rot_mat(:, :, i)';
    vtx0 = crystal_vtx(crystal.face{trace.fid(1)}, :);

    for k = 1:face_cnt
        n0 = trace_n(k);
        n1 = trace_n(k + 1);
        fn = face_norm(trace.fid(k), :);
        if n0 * n1 > 0
            % Refract
            [t, r, q] = compute_transit_factor(r, fn, n0, n1);
        elseif abs(n0) > 1 + 1e-6
            % Reflect in crystal
            [t, ~, ~] = compute_transit_factor(r, fn, abs(n0), 1);
            if isnan(t)
                t = 1;
            else
                t = 1 - t;
            end
            r = opt.reflect(r, fn);
            q = 1;
        else
            % Reflect out crystal
            [t, ~, ~] = compute_transit_factor(r, fn, 1, crystal.n);
            t = 1 - t;
            q = 1;
            r = opt.reflect(r, fn);
        end
        transit_factor(i) = transit_factor(i) * t * q;

        if k < length(trace.fid) && geo_factor(i) > 1e-8
            vtx1 = crystal_vtx(crystal.face{trace.fid(k + 1)}, :);
            [a, tmp_vtx] = face_intersection_factor(vtx0, vtx1, r);
            geo_factor(i) = geo_factor(i) * a;
            vtx0 = tmp_vtx;
        end
    end
end
transit_factor(isnan(transit_factor)) = 0;
geo_factor(isnan(geo_factor)) = 0;
end

% ================================================================================
function [t, ray_out, q_factor] = compute_transit_factor(ray_in, face_normal, n0, n1)
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

% ================================================================================
function [area_factor, vtx] = face_intersection_factor(vtx_entry, vtx_exit, ray_xyz)
% INPUT
%   vtx_entry:  entry face vertex
%   vtx_exit:   exit face vertex. it must be *CONVEX*
%   ray_xyz:    projection ray

e1 = vtx_exit(2, :) - vtx_exit(1, :);
e2 = vtx_exit(4, :) - vtx_exit(1, :);
e3 = cross(e1, e2);

uv0 = bsxfun(@minus, vtx_exit, vtx_exit(1, :)) / [e1; e2; e3];
uv0 = uv0(:, 1:2);
uvt = bsxfun(@minus, vtx_entry, vtx_exit(1, :)) / [e1; e2; -ray_xyz];
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
vtx = uv * [vtx_exit(2, :) - vtx_exit(1, :); vtx_exit(4, :) - vtx_exit(1, :)];
vtx = bsxfun(@plus, vtx, vtx_exit(1, :));

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
