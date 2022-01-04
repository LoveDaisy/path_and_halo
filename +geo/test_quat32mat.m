function test_quat32mat()
test_cases = {@suite1, @suite2, @suite3};
num = length(test_cases);

fprintf('Start testing for quat32mat...\n');
for i = 1:num
    fprintf('Testing suite %d/%d...\n', i, num);
    test_cases{i}();
    fprintf('suite %d/%d passed!\n', i, num);
end
end

% ================================================================================
function suite1()
quat3 = [.3, .5, - .7];
quat3 = quat3 / norm(quat3) * 0.7;
q = [sqrt(1 - norm(quat3)^2), quat3];

mat0 = quatrotate(q, eye(3))';

[mat, jac] = geo.quat32mat(quat3);
assert(all(abs(mat(:) - mat0(:)) < 1e-8));

h = 1e-4;
mat1 = geo.quat32mat(quat3 + [h, 0, 0]);
mat2 = geo.quat32mat(quat3 + [0, h, 0]);
mat3 = geo.quat32mat(quat3 + [0, 0, h]);
jac0 = cat(3, (mat1 - mat) / h, (mat2 - mat) / h, (mat3 - mat) / h);
assert(all(abs(jac(:) - jac0(:)) < 5e-4));
end

% ================================================================================
function suite2()
num = 100;
rng(1357);
q = randn(num, 4);
q(:, 1) = abs(q(:, 1)) + 0.01 * sqrt(sum(q.^2, 2));
q = bsxfun(@times, q, 1 ./ sqrt(sum(q.^2, 2)));
quat3 = q(:, 2:4);

mat0 = zeros(3, 3, num);
for i = 1:num
    mat0(:, :, i) = quatrotate(q(i, :), eye(3))';
end

[mat, jac] = geo.quat32mat(quat3);
assert(all(abs(mat(:) - mat0(:)) < 1e-8));

h = 1e-9;
mat1 = geo.quat32mat(bsxfun(@plus, quat3, [h, 0, 0]));
mat2 = geo.quat32mat(bsxfun(@plus, quat3, [0, h, 0]));
mat3 = geo.quat32mat(bsxfun(@plus, quat3, [0, 0, h]));
jac0 = cat(4, (mat1 - mat) / h, (mat2 - mat) / h, (mat3 - mat) / h);
jac0 = permute(jac0, [1, 2, 4, 3]);
assert(all(abs(jac(:) - jac0(:)) < 8e-4));
end

% ================================================================================
function suite3()
quat3 = [2, -1, 5];
quat3 = quat3 / norm(quat3);
q = [sqrt(1 - norm(quat3)^2), quat3];

mat0 = quatrotate(q, eye(3))';

[mat, jac] = geo.quat32mat(quat3);
assert(all(abs(mat(:) - mat0(:)) < 1e-8));

h = 1e-6;
mat1 = quatrotate([0, quat3 + [h, 0, 0]], eye(3))';
mat2 = quatrotate([0, quat3 + [0, h, 0]], eye(3))';
mat3 = quatrotate([0, quat3 + [0, 0, h]], eye(3))';
jac0 = cat(3, (mat1 - mat) / h, (mat2 - mat) / h, (mat3 - mat) / h);
assert(all(abs(jac(:) - jac0(:)) < 5e-4));
end
