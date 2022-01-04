function [mat, jac] = quat2mat(quat)
% Rotate matrix corresponding to quaternion [w, x, y, z].
% A rotation can be represented by a unit quaternion [w, xi, yj, zk],
% (see https://en.wikipedia.org/wiki/Rotation_matrix#Quaternion)
%
% INPUT
%   quat:       n*4, [w, x, y, z], a quaternion.
%
% OUTPUT
%   mat:        3*3*n, rotation matrix. It is for column vector, i.e. new_v = mat * v.
%   jac:        (3*3)*4*n, Jacobian, reshape first dimesion as matrix form

num = size(quat, 1);

[quat, jac_norm] = geo.normalize_vector(quat);

w = quat(:, 1);
x = quat(:, 2);
y = quat(:, 3);
z = quat(:, 4);

% The rotation matrix is (it is from MATLAB documants, seems different from wiki):
%   1 - 2*y^2 - 2*z^2,  2*x*y + 2*z*w,      2*x*z - 2*y*w
%   2*x*y - 2*z*w,      1 - 2*x^2 - 2*z^2,  2*y*z + 2*x*w
%   2*x*z + 2*y*w,      2*y*z - 2*x*w,      1 - 2*x^2 - 2*y^2

mat = [1 - 2 .* y.^2 - 2 .* z.^2, 2 .* x .* y - 2 .* z .* w, 2 .* x .* z + 2 .* y .* w, ...
    2 .* x .* y + 2 .* z .* w, 1 - 2 .* x.^2 - 2 .* z.^2, 2 .* y .* z - 2 .* x .* w, ...
        2 .* x .* z - 2 .* y .* w, 2 .* y .* z + 2 .* x .* w, 1 - 2 .* x.^2 - 2 .* y.^2];
mat = reshape(mat', [3, 3, num]);

jac = nan(3, 3, 4, num);
for i = 1:num
    jac1 = [0, 2 * z(i), -2 * y(i);
        -2 * z(i), 0, 2 * x(i);
        2 * y(i), -2 * x(i), 0];
    jac2 = [0, 2 * y(i), 2 * z(i);
        2 * y(i), -4 * x(i), 2 * w(i);
        2 * z(i), -2 * w(i), -4 * x(i)];
    jac3 = [-4 * y(i), 2 * x(i), -2 * w(i);
        2 * x(i), 0, 2 * z(i);
        2 * w(i), 2 * z(i), -4 * y(i)];
    jac4 = [-4 * z(i), 2 * w(i), 2 * x(i);
        -2 * w(i), -4 * z(i), 2 * y(i);
        2 * x(i), 2 * y(i), 0];

    temp_jac = [jac1(:), jac2(:), jac3(:), jac4(:)] * jac_norm(:, :, i);

    jac(:, :, 1, i) = reshape(temp_jac(:, 1), [3, 3]);
    jac(:, :, 2, i) = reshape(temp_jac(:, 2), [3, 3]);
    jac(:, :, 3, i) = reshape(temp_jac(:, 3), [3, 3]);
    jac(:, :, 4, i) = reshape(temp_jac(:, 4), [3, 3]);
end
end
