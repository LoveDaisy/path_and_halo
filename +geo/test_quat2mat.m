function test_quat2mat()
test_cases = {@suite1, @suite2};
num = length(test_cases);

fprintf('Start testing for quat2mat...\n');
for i = 1:num
    fprintf('Testing suite %d/%d...\n', i, num);
    test_cases{i}();
    fprintf('suite %d/%d passed!\n', i, num);
end
end

% ================================================================================
function suite1()
q3 = [.3, .5, - .7];
q3 = q3 / norm(q3) * 0.7;
q = [sqrt(1 - norm(q3)^2), q3];

mat0 = quatrotate(q, eye(3))';

[mat, jac] = geo.quat2mat(q);
assert(all(abs(mat(:) - mat0(:)) < 1e-12));

h = 1e-4;
mat1 = quatrotate(q + [h, 0, 0, 0], eye(3))';
mat2 = quatrotate(q + [0, h, 0, 0], eye(3))';
mat3 = quatrotate(q + [0, 0, h, 0], eye(3))';
mat4 = quatrotate(q + [0, 0, 0, h], eye(3))';
jac0 = cat(3, (mat1 - mat) / h, (mat2 - mat) / h, (mat3 - mat) / h, (mat4 - mat) / h);
assert(all(abs(jac(:) - jac0(:)) < 3e-4));
end

% ================================================================================
function suite2()
num = 100;
rng(1357);
q = randn(num, 4);

mat0 = zeros(3, 3, num);
for i = 1:num
    mat0(:, :, i) = quatrotate(q(i, :), eye(3))';
end

[mat, jac] = geo.quat2mat(q);
assert(all(abs(mat(:) - mat0(:)) < 1e-8));

h = 1e-9;
mat1 = geo.quat2mat(bsxfun(@plus, q, [h, 0, 0, 0]));
mat2 = geo.quat2mat(bsxfun(@plus, q, [0, h, 0, 0]));
mat3 = geo.quat2mat(bsxfun(@plus, q, [0, 0, h, 0]));
mat4 = geo.quat2mat(bsxfun(@plus, q, [0, 0, 0, h]));
jac0 = cat(4, (mat1 - mat) / h, (mat2 - mat) / h, (mat3 - mat) / h, (mat4 - mat) / h);
jac0 = permute(jac0, [1, 2, 4, 3]);
assert(all(abs(jac(:) - jac0(:)) < 8e-4));
end
