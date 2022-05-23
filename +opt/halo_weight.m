function [weight, status] = halo_weight(target_ll, fdf, config, axis_pdf, contour_info)
% High level function. Find a rotation contour for a given output LL

weight = 0;
status.fun_eval_cnt = 0;

[contour_store, inner_status] = opt.find_all_contours(target_ll, fdf, config, contour_info);
status.fun_eval_cnt = inner_status.fun_eval_cnt;

if isempty(contour_store)
    weight = -1;
    return;
end

for i = 1:length(contour_store)
    rot_contour = contour_store{i};
    % Compute weight
    [curr_w, curr_cmp, curr_rot] = opt.compute_contour_weight(rot_contour, axis_pdf, config);
    weight = weight + curr_w;
    contour_info.add_contour(rot_contour, curr_cmp, curr_rot);
end

end
