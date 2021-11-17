function [w, interp_s, interp_p, interp_rot] = ...
    compute_axis_rot_weight(curr_x, curr_jacob, axis_pdf, crystal, ...
    trace, sun_ll)

ray_in_xyz = ll2xyz_with_gradient(sun_ll);

p = axis_pdf(curr_x);
curr_det_j = zeros(size(curr_x, 1), 1);
face_area_factor = zeros(length(p), 1);
t_factor = zeros(length(p), 1);

for i = 1:length(p)
    curr_det_j(i) = det(curr_jacob(:, :, i) * curr_jacob(:, :, i)');
    tmp_rot_mat = rotz(90 + curr_x(i, 1)) * rotx(90 - curr_x(i, 2)) * rotz(curr_x(i, 3));
    
    tmp_crystal = crystal;
    tmp_crystal.face_norm = crystal.face_norm * tmp_rot_mat';
    tmp_crystal.vtx = crystal.vtx * tmp_rot_mat';
   
    tmp_face_area = entry_face_factor(tmp_crystal, ray_in_xyz);
    face_area_factor(i) = tmp_face_area(trace.fid(1));

    t_factor(i) = 1;
    trace_n = [1; trace.n];
    vtx0 = tmp_crystal.vtx(tmp_crystal.face{trace.fid(1)}, :);
    r = ray_in_xyz;
    for j = 1:length(trace.fid)
        [t, r] = transit_factor(r, ...
            tmp_crystal.face_norm(trace.fid(j), :), trace_n(j), trace_n(j+1));
        t_factor(i) = t_factor(i) * t;
        
        if j < length(trace.fid)
            vtx1 = tmp_crystal.vtx(tmp_crystal.face{trace.fid(j+1)}, :);
            [a, tmp_vtx] = face_intersection_factor(vtx0, vtx1, r);
            face_area_factor(i) = face_area_factor(i) * a;
            vtx0 = tmp_vtx;
        end
    end
end

[interp_rot, s, interp_s] = spline_interp_rot(curr_x);
tmp_det_j = exp(interp1(s, log(curr_det_j), interp_s, 'spline'));
tmp_face_factor = interp1(s, face_area_factor, interp_s, 'linear', 'extrap');
tmp_t_factor = exp(interp1(s, log(t_factor), interp_s, 'spline'));
interp_p = axis_pdf(interp_rot);
interp_p = [interp_p ./ tmp_det_j .* tmp_face_factor .* tmp_t_factor, ...
    interp_p, 1 ./ tmp_det_j, tmp_face_factor, tmp_t_factor];
w = sum((interp_p(1:end-1, 1) + interp_p(2:end, 1)) / 2 .* diff(interp_s));
end


function a = entry_face_factor(crystal, ray_in_xyz)
face_area = crystal.face_area .* (-crystal.face_norm * ray_in_xyz');
a = max(face_area, 0) / sum(max(face_area, 0));
end


function [t, ray_out] = transit_factor(ray_in, face_normal, n0, n1)
ray_out = refract_with_gradient(ray_in, face_normal, n0, n1);
cos_qi = -dot(ray_in, face_normal);
cos_qt = -dot(ray_out, face_normal);
Rs1 = abs((n0 * cos_qi - n1 * cos_qt) / (n0 * cos_qi + n1 * cos_qt))^2;
Rp1 = abs((n0 * cos_qt - n1 * cos_qi) / (n0 * cos_qt + n1 * cos_qi))^2;
t = (1 - (Rs1 + Rp1) / 2) * (cos_qi / cos_qt);
end


function [area_factor, vtx] = face_intersection_factor(vtx_q, vtx_0, ray_xyz)
% INPUT
%   vtx_q:      entry face vertex
%   vtx_0:      exit face vertex. it must be *CONVEX*
%   ray_xyz:    projection ray

uv0 = [0, 0; 1, 0; 1, 1; 0, 1];
uvt = bsxfun(@minus, vtx_q, vtx_0(1, :)) / ...
    [vtx_0(2, :) - vtx_0(1, :);
    vtx_0(4, :) - vtx_0(1, :);
    -ray_xyz];
uvq = uvt(:, 1:2);

mask_vtx_num = size(uv0, 1);
q_vtx_num = size(uvq, 1);
uv = nan(q_vtx_num * 2, 2);
for i1 = 1:mask_vtx_num
    i2 = mod(i1, mask_vtx_num) + 1;
    i3 = mod(i2, mask_vtx_num) + 1;

    z0 = det([uv0(i2, :) - uv0(i1, :); uv0(i3, :) - uv0(i2, :)]);
    uv_i = 1;
    for j1 = 1:q_vtx_num
        j2 = mod(j1, q_vtx_num) + 1;
        z1 = det([uv0(i2, :) - uv0(i1, :); uvq(j1, :) - uv0(i1, :)]);
        z2 = det([uv0(i2, :) - uv0(i1, :); uvq(j2, :) - uv0(i1, :)]);
        
        if z0 * z1 > 0 && z0 * z2 < 0
            % From inner to outer: record tow points
            uv(uv_i, :) = uvq(j1, :);
            ab = -(uv0(i1, :) - uvq(j1, :)) / [uv0(i2, :) - uv0(i1, :); uvq(j2, :) - uvq(j1, :)];
            uv(uv_i + 1, :) = uv0(i1, :) + ab(1) * (uv0(i2, :) - uv0(i1, :));
            uv_i = uv_i + 2;
        elseif z0 * z1 < 0 && z0 * z2 > 0
            % From outer to inner: record one point
            ab = -(uv0(i1, :) - uvq(j1, :)) / [uv0(i2, :) - uv0(i1, :); uvq(j2, :) - uvq(j1, :)];
            uv(uv_i, :) = uv0(i1, :) + ab(1) * (uv0(i2, :) - uv0(i1, :));
            uv_i = uv_i + 1;
        elseif z0 * z1 > 0 && z0 * z2 > 0
            % Two inner: record one point
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
    area_factor = polyarea(uv(:, 1), uv(:, 2));
else
    area_factor = 0;
end
end


function [interp_rot, s, interp_s] = spline_interp_rot(curr_x)
ds = sqrt(sum(diff(curr_x).^2, 2));
s = [0; cumsum(ds)];
num = size(curr_x, 1);

is_period = norm(curr_x(1, :) - curr_x(end, :)) < 1e-8;

ss = 0.05;
s_diff = inf;
iter_num = 1;
while s_diff > 1e-3 && iter_num < 3
    interp_s = [(s(1):ss:s(end))'; s];
    [c, ~, ic] = unique(interp_s);
    [interp_s, ics] = sort(c);
    [~, isc] = sort(ics);
    if is_period
        interp_rot = interp1([s; s(end) + s(2:end)], [curr_x; curr_x(2:end, :)], interp_s, 'spline');
        interp_rot = interp_rot(1:length(interp_s), :);
    else
        interp_rot = interp1(s, curr_x, interp_s, 'spline');
    end

    new_ds = sqrt(sum(diff(interp_rot).^2, 2));
    new_s = [0; cumsum(new_ds)];
    s_diff = abs(new_s(end) - s(end)) / s(end);
    
    s = new_s(isc(ic(end-num+1:end)));
    iter_num = iter_num + 1;
end
end
