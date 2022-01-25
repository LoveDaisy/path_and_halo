function [weight, status] = halo_weight(target_ll, fdf, config, axis_pdf, contour_info)
% High level function. Find a rotation contour for a given output LL

weight = -1;
status.fun_eval_cnt = 0;

% Find seed rotation
cand_rot = opt.find_cand_rot(config, target_ll, 'quat');
if isempty(cand_rot)
    return;
end

contour_h = 0.05;
reduce_eps = 0.05;
contour_info.add_candidate_seeds(cand_rot);

ray_out_xyz = geo.ll2xyz(target_ll);
weight = 0;
fun_eval_cnt = 0;
while ~isempty(cand_rot)
    % Find initial solution
    [init_rot, sol_status] = ode.find_solution(fdf, cand_rot(1, :), [ray_out_xyz, 1], 'eps', 1e-8);
    fun_eval_cnt = fun_eval_cnt + sol_status.fun_eval_cnt;
    cand_rot = cand_rot(2:end, :);
    if ~sol_status.finish || init_rot(1) < 0
        continue;
    end

    % Check if it locates on previous contour
    duplicated = false;
    for i = 1:contour_info.total_cnt
        [init_rot, dup_status] = geo.reduce_pts_polyline(contour_info.contour_store{i}, init_rot, ...
            'eps', reduce_eps);
        fun_eval_cnt = fun_eval_cnt + dup_status.fun_eval_cnt;
        if isempty(init_rot)
            duplicated = true;
            break;
        end
    end
    if duplicated
        continue;
    end

    % Find contour
    [rot_contour, contour_status] = ode.find_contour(fdf, init_rot, 'h', contour_h);
    fun_eval_cnt = fun_eval_cnt + contour_status.fun_eval_cnt;

    if isempty(rot_contour)
        continue;
    end

    % Reduce seeds
    [cand_rot, reduce_status] = geo.reduce_pts_polyline(rot_contour, cand_rot, ...
    'eps', config.dr * 2.5, 'jac_fun', fdf);
    fun_eval_cnt = fun_eval_cnt + reduce_status.fun_eval_cnt;

    % Compute weight
    [curr_w, curr_cmp, curr_rot] = opt.compute_contour_weight(rot_contour, axis_pdf, config);
    weight = weight + curr_w;
    contour_info.add_contour(rot_contour, curr_cmp, curr_rot);
end
status.fun_eval_cnt = fun_eval_cnt;
end
