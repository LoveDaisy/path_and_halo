function test_check_raypath()
test_cases = {@suite1};
num = length(test_cases);

fprintf('Start testing for check_raypath...\n');
for i = 1:num
    fprintf('Testing suite %d/%d...\n', i, num);
    test_cases{i}();
    fprintf('suite %d/%d passed!\n', i, num);
end
end

% ================================================================================
function suite1()
crystal = opt.make_prism_crystal(1);

valid_path = {[3, 5], [7, 5], [4, 2, 6], [3, 1, 6, 2], ...
    [1, 2, 1, 2, 1], [1, 3, 2, 5, 1], [1, 2, 3, 5, 1], [1, 2, 3, 4, 1]};
for i = 1:length(valid_path)
    fprintf(' valid path %d/%d... ', i, length(valid_path));
    assert(opt.check_raypath(crystal, valid_path{i}));
    fprintf('pass!\n');
end

invalid_path = {[3, 4], [1, 3, 4], [5, 6, 3], [1, 3, 1, 3], [1, 2, 3, 6], [3, 1, 2, 6], ...
    [3, 1, 3, 1, 3]};
for i = 1:length(invalid_path)
    fprintf(' invalid path %d/%d... ', i, length(invalid_path));
    assert(~opt.check_raypath(crystal, invalid_path{i}));
    fprintf('pass!\n');
end
end