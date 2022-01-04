function [mat, jac] = quat32mat(quat3)
% Rotate matrix corresponding to quat3-representation, i.e. [x, y, z] for imaginary part of a unit quaternion.
% A rotation can be represented by a unit quaternion [w, xi, yj, zk],
% (see https://en.wikipedia.org/wiki/Rotation_matrix#Quaternion)
%
% INPUT
%   quat3:      n*3, [x, y, z], the imaginary party of a unit quaternion
%
% OUTPUT
%   mat:        3*3*n, rotation matrix. It is for column vector, i.e. new_v = mat * v.
%   jac:        (3*3)*3*n, Jacobian, reshape first dimesion as matrix form

p = inputParser;
p.addRequired('quat3', @(x) validateattributes(x, {'numeric'}, {'2d', 'ncols', 3}));
p.parse(quat3);

num = size(quat3, 1);

n2 = sum(quat3.^2, 2);
invalid_idx = n2 >= 1;

jac_norm = nan(3, 3, num);
if any(invalid_idx)
    [quat3(invalid_idx, :), tmp_jac_norm] = geo.normalize_vector(quat3(invalid_idx, :));
    jac_norm(:, :, invalid_idx) = tmp_jac_norm;
end
for i = 1:num
    if ~invalid_idx(i)
        jac_norm(:, :, i) = eye(3);
    end
end

w = sqrt(max(1 - n2, 0));
x = quat3(:, 1);
y = quat3(:, 2);
z = quat3(:, 3);

% The rotation matrix is (it is from MATLAB documants, seems different from wiki):
%   1 - 2*y^2 - 2*z^2,  2*x*y + 2*z*w,      2*x*z - 2*y*w
%   2*x*y - 2*z*w,      1 - 2*x^2 - 2*z^2,  2*y*z + 2*x*w
%   2*x*z + 2*y*w,      2*y*z - 2*x*w,      1 - 2*x^2 - 2*y^2

w1 = -x ./ w;
w1(invalid_idx) = 0;
w2 = -y ./ w;
w2(invalid_idx) = 0;
w3 = -z ./ w;
w3(invalid_idx) = 0;

mat = nan(3, 3, num);
jac = nan(3, 3, 3, num);
for i = 1:num
    m = [1 - 2 * y(i)^2 - 2 * z(i)^2, 2 * x(i) * y(i) + 2 * z(i) * w(i), 2 * x(i) * z(i) - 2 * y(i) * w(i);
        2 * x(i) * y(i) - 2 * z(i) * w(i), 1 - 2 * x(i)^2 - 2 * z(i)^2, 2 * y(i) * z(i) + 2 * x(i) * w(i);
        2 * x(i) * z(i) + 2 * y(i) * w(i), 2 * y(i) * z(i) - 2 * x(i) * w(i), 1 - 2 * x(i)^2 - 2 * y(i)^2];
    jac1 = [0, 2 * y(i) + 2 * z(i) * w1(i), 2 * z(i) - 2 * y(i) * w1(i);
        2 * y(i) - 2 * z(i) * w1(i), -4 * x(i), 2 * w(i) + 2 * x(i) * w1(i);
        2 * z(i) + 2 * y(i) * w1(i), -2 * w(i) - 2 * x(i) * w1(i), -4 * x(i)];
    jac2 = [-4 * y(i), 2 * x(i) + 2 * z(i) * w2(i), -2 * w(i) - 2 * y(i) * w2(i);
        2 * x(i) - 2 * z(i) * w2(i), 0, 2 * z(i) + 2 * x(i) * w2(i);
        2 * w(i) + 2 * y(i) * w2(i), 2 * z(i) - 2 * x(i) * w2(i), -4 * y(i)];
    jac3 = [-4 * z(i), 2 * w(i) + 2 * z(i) * w3(i), 2 * x(i) - 2 * y(i) * w3(i);
        -2 * w(i) - 2 * z(i) * w3(i), -4 * z(i), 2 * y(i) + 2 * x(i) * w3(i);
        2 * x(i) + 2 * y(i) * w3(i), 2 * y(i) - 2 * x(i) * w3(i), 0];

    mat(:, :, i) = m;
    if invalid_idx(i)
        tmp_jac = [jac1(:), jac2(:), jac3(:)] * jac_norm(:, :, i);
        jac(:, :, 1, i) = reshape(tmp_jac(:, 1), [3, 3]);
        jac(:, :, 2, i) = reshape(tmp_jac(:, 2), [3, 3]);
        jac(:, :, 3, i) = reshape(tmp_jac(:, 3), [3, 3]);
    else
        jac(:, :, 1, i) = jac1;
        jac(:, :, 2, i) = jac2;
        jac(:, :, 3, i) = jac3;
    end
end
end
