function [q, jac] = llr2quat(llr)
% Convert a 3D rotation from [longitude, latitude, roll] to quaternion [w, x, y, z]
% A rotation can be represented by an axis direction and a roll around the axis.
% Specially in this project, we use a direction of crystal main axis and a roll around it, which
% is a little different from nomal axis-angle representation
% (see https://en.wikipedia.org/wiki/Axis%E2%80%93angle_representation)
% The rotation matrix is (for column vector, rotation matrix is left-multiplied):
%   rotz(90 + longitude) * rotx(90 - latitude) * rotz(roll)
% which means:
%   1. apply roll around z-axis to object,
%   2. apply 90 - latitude around x-axis to object,
%   3. apply 90 + longitude around z-axis to object.
%
% INPUT
%   llr:        n*3, [longitude, latitude, roll], in degree
%
% OUTPUT
%   q:          n*4, quaternion. Apply quaternion using quatrotate(q, v) to row vector has the same
%               result with using matrix like v * mat' (v is row vector thus the matrix is transposed
%               and right-multiplied)
%   jac :       4*3*n, Jacobian of this transformation

p = inputParser;
p.addRequired('llr', @(x) validateattributes(x, {'numeric'}, {'2d', 'ncols', 3}));
p.parse(llr);

% Construct three quaternion corresponding to three rotation matrix
q1 = [cosd(llr(:, 3) / 2), -sind(llr(:, 3) / 2) * [0, 0, 1]];
q2 = [cosd(45 - llr(:, 2) / 2), -sind(45 - llr(:, 2) / 2) * [1, 0, 0]];
q3 = [cosd(45 + llr(:, 1) / 2), -sind(45 + llr(:, 1) / 2) * [0, 0, 1]];
qm = quatmultiply(q1, q2);
q = quatmultiply(qm, q3);
[q, jac_norm] = geo.normalize_vector(q);

num = size(llr, 1);
jac = nan(4, 3, num);
for i = 1:num
    j1 = [0, 0, -sind(llr(i, 3) / 2);
        0, 0, 0;
        0, 0, 0;
        0, 0, -cosd(llr(i, 3) / 2)] * pi / 180/2;
    j2 = [0, sind(45 - llr(i, 2) / 2), 0;
        0, cosd(45 - llr(i, 2) / 2), 0;
        0, 0, 0;
        0, 0, 0] * pi / 180/2;
    j3 = [-sind(45 + llr(i, 1) / 2), 0, 0;
        0, 0, 0;
        0, 0, 0;
        -cosd(45 + llr(i, 1) / 2), 0, 0] * pi / 180/2;
    jm = [q2(i, 1), -q2(i, 2), -q2(i, 3), -q2(i, 4);
        q2(i, 2), q2(i, 1), q2(i, 4), -q2(i, 3);
        q2(i, 3), -q2(i, 4), q2(i, 1), q2(i, 2);
        q2(i, 4), q2(i, 3), -q2(i, 2), q2(i, 1)] * j1 + ...
        [q1(i, 1), -q1(i, 2), -q1(i, 3), -q1(i, 4);
        q1(i, 2), q1(i, 1), -q1(i, 4), q1(i, 3);
        q1(i, 3), q1(i, 4), q1(i, 1), -q1(i, 2);
        q1(i, 4), -q1(i, 3), q1(i, 2), q1(i, 1)] * j2;
    jq = [q3(i, 1), -q3(i, 2), -q3(i, 3), -q3(i, 4);
        q3(i, 2), q3(i, 1), q3(i, 4), -q3(i, 3);
        q3(i, 3), -q3(i, 4), q3(i, 1), q3(i, 2);
        q3(i, 4), q3(i, 3), -q3(i, 2), q3(i, 1)] * jm + ...
        [qm(i, 1), -qm(i, 2), -qm(i, 3), -qm(i, 4);
        qm(i, 2), qm(i, 1), -qm(i, 4), qm(i, 3);
        qm(i, 3), qm(i, 4), qm(i, 1), -qm(i, 2);
        qm(i, 4), -qm(i, 3), qm(i, 2), qm(i, 1)] * j3;
    jac(:, :, i) = jac_norm(:, :, i) * jq;
end
end
