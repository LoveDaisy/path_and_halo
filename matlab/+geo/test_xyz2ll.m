function test_xyz2ll()
test_cases = {@suite1, @suite2};
num = length(test_cases);

fprintf('Start testing for xyz2ll...\n');
for i = 1:num
    fprintf('Testing suite %d/%d...\n', i, num);
    test_cases{i}();
    fprintf('suite %d/%d passed!\n', i, num);
end
end

% ================================================================================
function suite1()
xyz = [1, 2, 3];
h = 1e-6;
ll0 = [atan2d(xyz(2), xyz(1)), asind(xyz(3) / norm(xyz))];
tmp_ll1 = [atan2d(xyz(2), xyz(1) + h), asind(xyz(3) / norm(xyz + [h, 0, 0]))];
tmp_ll2 = [atan2d(xyz(2) + h, xyz(1)), asind(xyz(3) / norm(xyz + [0, h, 0]))];
tmp_ll3 = [atan2d(xyz(2), xyz(1)), asind((xyz(3) + h) / norm(xyz + [0, 0, h]))];
jac0 = [tmp_ll1 - ll0; tmp_ll2 - ll0; tmp_ll3 - ll0]' / h;

[ll, jac] = geo.xyz2ll(xyz);
assert(all(abs(ll(:) - ll0(:)) < 1e-6));
assert(all(abs(jac(:) - jac0(:)) < 1e-5));
end

% ================================================================================
function suite2()
num = 100;
rng(1234);
xyz = randn(num, 3);
h = 1e-8;
ll0 = [atan2d(xyz(:, 2), xyz(:, 1)), asind(xyz(:, 3) ./ sqrt(sum(xyz.^2, 2)))];
jac0 = zeros(2, 3, num);
for i = 1:num
    tmp_ll1 = [atan2d(xyz(i, 2), xyz(i, 1) + h), asind(xyz(i, 3) / norm(xyz(i, :) + [h, 0, 0]))];
    tmp_ll2 = [atan2d(xyz(i, 2) + h, xyz(i, 1)), asind(xyz(i, 3) / norm(xyz(i, :) + [0, h, 0]))];
    tmp_ll3 = [atan2d(xyz(i, 2), xyz(i, 1)), asind((xyz(i, 3) + h) / norm(xyz(i, :) + [0, 0, h]))];
    jac0(:, :, i) = [tmp_ll1 - ll0(i, :); tmp_ll2 - ll0(i, :); tmp_ll3 - ll0(i, :)]' / h;
end

[ll, jac] = geo.xyz2ll(xyz);
assert(all(abs(ll(:) - ll0(:)) < 1e-6));
assert(all(abs(jac(:) - jac0(:)) < 5e-4));
end
