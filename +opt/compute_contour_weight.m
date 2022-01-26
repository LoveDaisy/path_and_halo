function [weight, cmp_interp, llr_interp] = compute_contour_weight(rot_contour, axis_pdf, config)
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
%   cmp_interp:         m*6, [s, w, axis_p, det_jac, geo_factor, transit_factor]
%   llr_interp:         m*3, interpolated rotation, in LLR space

dim = size(rot_contour, 2);
pdf_th = 1e-10;

cmp0 = compute_weight_components(rot_contour, axis_pdf, config);

% Calculate length of polyline as initial value
len = sum(sqrt(sum(diff(rot_contour).^2, 2)));
interp_step = len * 0.01;

% This while loop makes sure there are enough points with high PDF
while true
    % Interpolate rotations
    [rot_interp, s_interp, s0, s0_idx] = geo.interp_curve(rot_contour, interp_step);
    interp_num = size(rot_interp, 1);

    % Convert to LLR space
    if dim == 4
        llr_interp = geo.quat2llr(rot_interp);
        diff_s = sqrt(sum(diff(llr_interp).^2, 2));
        s_interp = [0; cumsum(diff_s)];
        discontinuity_idx = find(diff_s > 150) + 1;
        for i = length(discontinuity_idx):-1:1
            idx = discontinuity_idx(i);
            d1 = diff_s(idx - 1);
            if idx > 2
                d0 = diff_s(idx - 2);
            else
                d0 = 0;
            end
            s_interp(idx:end) = s_interp(idx:end) - d1 + d0;
        end
        s0 = s_interp(s0_idx);
    elseif dim == 3
        llr_interp = rot_interp;
    else
        error('Input rotation must have dimesion of 3 or 4!');
    end

    interp_pdf = axis_pdf(llr_interp) .* cosd(llr_interp(:, 2));
    if max(interp_pdf) > pdf_th && sum(interp_pdf >= pdf_th) < 20 && interp_step > len * 0.001
        interp_step = interp_step * 0.5;
    else
        break;
    end
end

% Interpolate components as initial value
cmp_interp = nan(length(s_interp), 6);
cmp_interp(:, 1) = s_interp;
cmp_interp(:, 3) = interp_pdf;
cmp_interp(:, 4) = exp(interp1(s0, log(cmp0(:, 2)), s_interp, 'linear', 'extrap'));
cmp_interp(:, 5) = interp1(s0, cmp0(:, 3), s_interp, 'linear', 'extrap');
cmp_interp(:, 6) = exp(interp1(s0, log(cmp0(:, 4)), s_interp, 'linear', 'extrap'));

% Find out those with large probability && (small geometry factor || small transit_factor)
interest_idx = interp_pdf >= pdf_th & (cmp_interp(:, 5) <= 1e-1 | cmp_interp(:, 6) <= 1e-1);
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
    cmp_interp(i1:i2, 3:6) = compute_weight_components(llr_interp(i1:i2, :), axis_pdf, config);
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
det_jac = jacobian_factor(llr, config);
entry_factor = entry_face_factor(llr, config);
[transit_factor, geo_factor] = transit_geo_factor(llr, config);

cmp = [axis_prob, 1 ./ det_jac, entry_factor .* geo_factor, transit_factor];
end

% ================================================================================
function det_jac = jacobian_factor(llr, config)
% Compute deteminant of Jacobian.
% What does this Jacobian mean? We simply imagine this, when a point moves within a small
% region orthogonal to contour (for LLR representation, the small region is a 2D subspace),
% the output will also moves with in a small region.
% And clearly, the determinant of Jacobian, is the ratio of these two reagion volumn.

rot_num = size(llr, 1);
ray_in_ll = [config.sun_ll(1) + 180, -config.sun_ll(2)];
fdf = @(rot) opt.crystal_system(rot, ray_in_ll, config.crystal, config.trace);

[~, jac] = fdf(llr);

min_det = 1e-8;
jac_rank_tol = 1e-10;

% Let [u, s, v] = svd(jac), i.e. u * s * v' = jac.
% If rank(jac) = k, then there are only k values greater than 0 in s.
% We devide u, s, v into blocks,
% u = [u_o, u_n],     v = [v_o, v_n]
%      m*k  m*(m-k)        n*k  n*(n-k)
% Then jac = u_o * s(1:k, 1:k) * v_o'
% We will find that v_n is null space of jac: jac * v_n = u_o * s(1:k, 1:k) * v_o' * v_n = 0,
% and thus v_n is the tangent direction of contour.
% Let dx be in the orthogonal complement to v_n, i.e. perpendicular to contour. It can be
% expressed as dx = v_o * alpha. So,
%   |dy|^2 = dx' * jac' * jac * dx
%          = alpha' * v_o' * v_o * s_o * u_o' * u_o * s_o * v_o' * v_o * alpha
%          = alpha' * s_o^2 * alpha
% So the area factor is det(s_o)
det_jac = nan(size(llr, 1), 1);
for i = 1:rot_num
    curr_jac = jac(:, :, i);
    s = svd(curr_jac);
    jac_rank = sum(s > jac_rank_tol);
    det_jac(i) = max(prod(s(1:jac_rank)), min_det);
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
            [a, tmp_vtx] = geo.project_polygon_intersection(vtx0, vtx1, r);
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
