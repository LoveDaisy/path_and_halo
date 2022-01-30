function test_init_config()
test_cases = {@suite1};
num = length(test_cases);

fprintf('Start testing for init_config...\n');
for i = 1:num
    fprintf('Testing suite %d/%d...\n', i, num);
    test_cases{i}();
    fprintf('suite %d/%d passed!\n', i, num);
end
end

% ================================================================================
function suite1()
crystal = opt.make_prism_crystal(1);
trace.fid = [3; 5];
sun_ll = [0, 25];

config = opt.init_config(crystal, trace, sun_ll, 'GridLevel', 2);

hist_n0 = [8, 130, 118, 98, 68, 72, 72, 54, 28, 42, 44, 22];
hist_edge0 = [20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44];
[n, edge] = histcounts(config.bending_angle);

assert(all(abs(hist_n0(:) - n(:)) < 1e-8));
assert(all(abs(hist_edge0(:) - edge(:)) < 1e-8));
end
