function test_find_seed_rot()
test_cases = {@suite1};
num = length(test_cases);

fprintf('Start testing for find_seed_rot...\n');
for i = 1:num
    fprintf('Testing suite %d/%d...\n', i, num);
    test_cases{i}();
    fprintf('suite %d/%d passed!\n', i, num);
end
end

% ================================================================================
function suite1()
sun_ll = [0, 10];
target_ll = [0.5, -12];

crystal = opt.make_prism_crystal(1);
trace.fid = [1; 3; 2; 4; 5; 1];

config_cache_file = 'test_config_1_132451_0+10_3.mat';
if exist(config_cache_file, 'file')
    load(config_cache_file);
else
    config = opt.init_config(crystal, trace, sun_ll, 3);
    save(config_cache_file, 'config');
end

ray_in_ll = [sun_ll(1) + 180, -sun_ll(2)];
fdf = @(rot) opt.crystal_system(rot, ray_in_ll, config.crystal, config.trace);

[seed_llr, status] = opt.find_seed_rot(config, target_ll, 'llr');
for i = 1:size(seed_llr, 1)
    assert(status.fun_eval_cnt < 300);

    tmp_ll = fdf(seed_llr(i, :));
    assert(norm(tmp_ll - target_ll) < 1e-8);
end

[seed_quat, status] = opt.find_seed_rot(config, target_ll, 'quat');
for i = 1:size(seed_quat, 1)
    assert(status.fun_eval_cnt < 300);

    tmp_ll = fdf(seed_quat(i, :));
    assert(norm(tmp_ll - [target_ll, 1]) < 1e-8);
end
end
