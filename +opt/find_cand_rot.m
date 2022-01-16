function cand_rot = find_cand_rot(config, target_ll, type)
out_ll_diff = acosd(geo.ll2xyz(config.out_ll) * geo.ll2xyz(target_ll)');
[out_ll_diff, sort_idx] = sort(out_ll_diff);

cand_idx = sort_idx(out_ll_diff < config.dr * 1.5);
cand_llr = config.axis_llr_store(cand_idx, :);
cand_quat = config.axis_quat_store(cand_idx, :);

if strcmpi(type, 'llr')
    cand_rot = cand_llr;
elseif strcmpi(type, 'quat')
    cand_rot = cand_quat;
else
    error('unknown rotation type!');
end
end
