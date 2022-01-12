function test_rotation_representation()
test_cases = {@suite1, @suite2, @suite3, @suite4};
num = length(test_cases);

fprintf('Start testing for conversion between rotation representations...\n');
for i = 1:num
    fprintf('Testing suite %d/%d...\n', i, num);
    test_cases{i}();
    fprintf('suite %d/%d passed!\n', i, num);
end
end

% ================================================================================
function suite1()
% Quaternion to matrix
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
% Quaternion to matrix
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

% ================================================================================
function suite3()
% LLR to quaternion
lon = 10;
lat = 35;
roll = 60;
mat = rotz(90 + lon) * rotx(90 - lat) * rotz(roll);

v = [3, 4, 1.2];
v_new0 = v * mat';

[q, jac] = geo.llr2quat([lon, lat, roll]);
v_new = quatrotate(q, v);

assert(all(abs(v_new - v_new0) < 1e-6));

d = 1e-4;
[q1, ~] = geo.llr2quat([lon, lat, roll] + [d, 0, 0]);
[q2, ~] = geo.llr2quat([lon, lat, roll] + [0, d, 0]);
[q3, ~] = geo.llr2quat([lon, lat, roll] + [0, 0, d]);

jac0 = [q1 - q; q2 - q; q3 - q]' / d;
assert(all(abs(jac0(:) - jac(:)) < 1e-8));
end

% ================================================================================
function suite4()
% LLR to matrix
lon = 10;
lat = 35;
roll = 60;
mat0 = rotz(90 + lon) * rotx(90 - lat) * rotz(roll);

[mat, jac] = geo.llr2mat([lon, lat, roll]);
assert(all(abs(mat(:) - mat0(:)) < 1e-10));

d = 1e-3;
mat1 = geo.llr2mat([lon, lat, roll] + [d, 0, 0]);
mat2 = geo.llr2mat([lon, lat, roll] + [0, d, 0]);
mat3 = geo.llr2mat([lon, lat, roll] + [0, 0, d]);
jac0 = cat(3, (mat1 - mat0) / d, (mat2 - mat0) / d, (mat3 - mat0) / d);
assert(all(abs(jac(:) - jac0(:)) < 5e-7));
end
