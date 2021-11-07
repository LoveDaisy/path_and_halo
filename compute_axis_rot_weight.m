function [w, interp_s, interp_p, interp_rot] = ...
    compute_axis_rot_weight(curr_x, curr_jacob, axis_pdf, all_face_norm, all_face_area, ...
    entry_face_idx, sun_ll)

ray_in_xyz = ll2xyz_with_gradient(sun_ll);

p = axis_pdf(curr_x);
curr_det_j = zeros(size(curr_x, 1), 1);
face_area_factor = zeros(length(p), 1);

for j = 1:length(p)
    curr_det_j(j) = det(curr_jacob(:, :, j) * curr_jacob(:, :, j)');
    tmp_rot_mat = rotz(90 + curr_x(j, 1)) * rotx(90 - curr_x(j, 2)) * rotz(curr_x(j, 3));
    tmp_face_norm = all_face_norm * tmp_rot_mat';
    tmp_face_area = all_face_area .* (-tmp_face_norm * ray_in_xyz');
    tmp_face_area = tmp_face_area(tmp_face_area > 0);
    face_area_factor(j) = max(tmp_face_area(entry_face_idx), 0) / sum(max(tmp_face_area, 0));
end
ds = sqrt(sum(diff(curr_x).^2, 2));
s = [0; cumsum(ds)];

interp_s = (s(1):0.05:s(end))';
interp_rot = interp1(s, curr_x, interp_s, 'spline');
tmp_det_j = interp1(s, curr_det_j, interp_s, 'spline');
tmp_face_factor = interp1(s, face_area_factor, interp_s, 'linear');
interp_p = axis_pdf(interp_rot);
interp_p = interp_p ./ tmp_det_j .* tmp_face_factor;
w = sum((interp_p(1:end-1) + interp_p(2:end)) / 2 .* diff(interp_s));
end