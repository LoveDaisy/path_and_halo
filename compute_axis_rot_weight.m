function [w, interp_s, interp_p, interp_rot] = ...
    compute_axis_rot_weight(curr_x, curr_jacob, axis_pdf, crystal, ...
    trace, sun_ll)

ray_in_xyz = ll2xyz_with_gradient(sun_ll);

p = axis_pdf(curr_x);
curr_det_j = zeros(size(curr_x, 1), 1);
face_area_factor = zeros(length(p), 1);

for j = 1:length(p)
    curr_det_j(j) = det(curr_jacob(:, :, j) * curr_jacob(:, :, j)');
    tmp_rot_mat = rotz(90 + curr_x(j, 1)) * rotx(90 - curr_x(j, 2)) * rotz(curr_x(j, 3));
    tmp_face_norm = crystal.face_norm * tmp_rot_mat';
    tmp_face_area = crystal.face_area .* (-tmp_face_norm * ray_in_xyz');
    tmp_face_area = tmp_face_area(tmp_face_area > 0);
    face_area_factor(j) = max(tmp_face_area(trace.fid(1)), 0) / sum(max(tmp_face_area, 0));
end

[interp_rot, s, interp_s] = spline_interp_rot(curr_x);
tmp_det_j = exp(interp1(s, log(curr_det_j), interp_s, 'spline'));
tmp_face_factor = interp1(s, face_area_factor, interp_s, 'linear', 'extrap');
interp_p = axis_pdf(interp_rot);
interp_p = [interp_p ./ tmp_det_j .* tmp_face_factor, interp_p, 1 ./ tmp_det_j, tmp_face_factor];
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
