function [w, interp_s, interp_p, interp_rot] = ...
    compute_axis_rot_weight(curr_x, curr_jacob, axis_pdf, crystal, ...
    trace, sun_ll)

ray_in_xyz = ll2xyz_with_gradient(sun_ll);

p = axis_pdf(curr_x);
curr_det_j = zeros(size(curr_x, 1), 1);
face_area_factor = zeros(length(p), 1);
t_factor = zeros(length(p), 1);

for j = 1:length(p)
    curr_det_j(j) = det(curr_jacob(:, :, j) * curr_jacob(:, :, j)');
    tmp_rot_mat = rotz(90 + curr_x(j, 1)) * rotx(90 - curr_x(j, 2)) * rotz(curr_x(j, 3));
    
    tmp_crystal = crystal;
    tmp_crystal.face_norm = crystal.face_norm * tmp_rot_mat';
    tmp_crystal.vtx = crystal.vtx * tmp_rot_mat';
    
    tmp_face_area = crystal.face_area .* (-tmp_crystal.face_norm * ray_in_xyz');
    tmp_face_area = tmp_face_area(tmp_face_area > 0);
    face_area_factor(j) = max(tmp_face_area(trace.fid(1)), 0) / sum(max(tmp_face_area, 0));

    tmp_face_norm = tmp_crystal.face_norm(trace.fid(1), :);
    tmp_refract_ray = refract_with_gradient(ray_in_xyz, tmp_face_norm, 1, trace.n(1));
    cos_qi = -dot(ray_in_xyz, tmp_face_norm);
    cos_qt = -dot(tmp_refract_ray, tmp_face_norm);
    Rs1 = abs((cos_qi - trace.n(1) * cos_qt) / (cos_qi + trace.n(1) * cos_qt))^2;
    Rp1 = abs((cos_qt - trace.n(1) * cos_qi) / (cos_qt + trace.n(1) * cos_qi))^2;
    T1 = (1 - (Rs1 + Rp1) / 2) * (cos_qi / cos_qt);
   
    tmp_face_norm = tmp_crystal.face_norm(trace.fid(end), :);
    tmp_exit_ray = refract_with_gradient(tmp_refract_ray, tmp_face_norm, trace.n(1), 1);
    cos_qi = dot(tmp_refract_ray, tmp_face_norm);
    cos_qt = dot(tmp_exit_ray, tmp_face_norm);
    Rs2 = abs((trace.n(1) * cos_qi - cos_qt) / (trace.n(1) * cos_qi + cos_qt))^2;
    Rp2 = abs((trace.n(1) * cos_qt - cos_qi) / (trace.n(1) * cos_qt + cos_qi))^2;
    T2 = (1 - (Rs2 + Rp2) / 2) * (cos_qi / cos_qt);

    t_factor(j) = T2 * T1;
    
    tmp_entry_vtx = tmp_crystal.vtx(tmp_crystal.face{trace.fid(1)}, :);
    tmp_exit_vtx = tmp_crystal.vtx(tmp_crystal.face{trace.fid(end)}, :);
    uvt = bsxfun(@minus, tmp_entry_vtx, tmp_exit_vtx(1, :)) / ...
        [tmp_exit_vtx(2, :) - tmp_exit_vtx(1, :);
        tmp_exit_vtx(4, :) - tmp_exit_vtx(1, :);
        -tmp_refract_ray];
    tmp_proj_vtx = tmp_entry_vtx + uvt(:, 3) * tmp_refract_ray;
    uv0 = [0, 0; 1, 0; 1, 1; 0, 1];
    [ux, vx] = polyxpoly([uvt(:, 1); uvt(1, 1)], [uvt(:, 2); uvt(1, 2)], ...
        [uv0(:, 1); uv0(1, 1)], [uv0(:, 2); uv0(1, 2)]);
    tmp_in0 = inpolygon(uvt(:, 1), uvt(:, 2), uv0(:, 1), uv0(:, 2));
    tmp_int = inpolygon(uv0(:, 1), uv0(:, 2), uvt(:, 1), uvt(:, 2));
    quv = unique([uvt(tmp_in0, 1:2); uv0(tmp_int, :); ux, vx], 'rows');
    if size(quv, 1) > 2
        tmp_idx = convhull(quv(:, 1), quv(:, 2));
        tmp_area = polyarea(quv(tmp_idx, 1), quv(tmp_idx, 2));
    else
        tmp_area = 0;
    end
    
    face_area_factor(j) = face_area_factor(j) * tmp_area;
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
