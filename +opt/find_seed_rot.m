function [seed_llr, seed_quat, status] = find_seed_rot(config, target_ll)
out_ll_diff = acosd(geo.ll2xyz(config.out_ll) * geo.ll2xyz(target_ll)');
cand_idx = out_ll_diff < config.dr;
cand_cnt = sum(cand_idx);
cand_llr = config.axis_llr_store(cand_idx, :);
cand_quat = config.axis_quat_store(cand_idx, :);

sun_ll = config.sun_ll;
ray_in_ll = [sun_ll(1) + 180, -sun_ll(2)];
fdf = @(rot) opt.crystal_system(rot, ray_in_ll, config.crystal, config.trace);
fun_eval_cnt = 0;

out_eps = 1e-8;
in_eps = 1e-4;
seed_llr = nan(size(cand_llr));
idx = 0;
for i = 1:cand_cnt
    [tmp_llr, tmp_status] = ode.find_solution_fdf(fdf, cand_llr(i, :), target_ll, 'eps', out_eps);
    fun_eval_cnt = fun_eval_cnt + tmp_status.fun_eval_cnt;
    tmp_seed_diff = sqrt(sum(bsxfun(@minus, seed_llr, tmp_llr).^2, 2));
    if tmp_status.finish && (min(tmp_seed_diff) > in_eps || all(isnan(tmp_seed_diff)))
        idx = idx + 1;
        seed_llr(idx, :) = tmp_llr;
    end
end
seed_llr = seed_llr(1:idx, :);

out_eps = 1e-8;
in_eps = 1e-6;
seed_quat = nan(size(cand_quat));
idx = 0;
for i = 1:cand_cnt
    [tmp_quat, tmp_status] = ode.find_solution_fdf(fdf, cand_quat(i, :), [target_ll, 1], 'eps', out_eps);
    fun_eval_cnt = fun_eval_cnt + tmp_status.fun_eval_cnt;
    tmp_seed_diff = sqrt(sum(bsxfun(@minus, seed_quat, tmp_quat).^2, 2));
    if tmp_status.finish && (min(tmp_seed_diff) > in_eps || all(isnan(tmp_seed_diff)))
        idx = idx + 1;
        seed_quat(idx, :) = tmp_quat;
    end
end
seed_quat = seed_quat(1:idx, :);

status.fun_eval_cnt = fun_eval_cnt;
end