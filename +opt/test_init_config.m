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

config = opt.init_config(crystal, trace, sun_ll, 2);
end
