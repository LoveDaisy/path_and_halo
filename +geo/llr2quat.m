function [q, grad_q] = llr2quat(llr)
% INPUT

% The rotation matrix is:
%  rotz(90 + llr(1)) * rotx(90 - llr(2)) * rotz(llr(3))
q1 = [cosd(llr(:, 3) / 2), -sind(llr(:, 3) / 2) * [0, 0, 1]];
q2 = [cosd(45 - llr(:, 2) / 2), -sind(45 - llr(:, 2) / 2) * [1, 0, 0]];
q3 = [cosd(45 + llr(:, 1) / 2), -sind(45 + llr(:, 1) / 2) * [0, 0, 1]];
qm = quatmultiply(q1, q2);
q = quatmultiply(qm, q3);
[q, grad_norm] = geo.normalize_vector(q);

num = size(llr, 1);
grad_q = nan(4, 3, num);
for i = 1:num
    j1 = [0, 0, -sind(llr(i, 3) / 2);
        0, 0, 0;
        0, 0, 0;
        0, 0, -cosd(llr(i, 3) / 2)] * pi / 180 / 2;
    j2 = [0, sind(45 - llr(i, 2) / 2), 0;
        0, cosd(45 - llr(i, 2) / 2), 0;
        0, 0, 0;
        0, 0, 0] * pi / 180 / 2;
    j3 = [-sind(45 + llr(i, 1) / 2), 0, 0;
        0, 0, 0;
        0, 0, 0;
        -cosd(45 + llr(i, 1) / 2), 0, 0] * pi / 180 / 2;
    jm = [q2(i, 1), -q2(i, 2), -q2(i, 3), -q2(i, 4);
        q2(i, 2), q2(i, 1), q2(i, 4), -q2(i, 3);
        q2(i, 3), -q2(i, 4), q2(i, 1), q2(i, 2);
        q2(i, 4), q2(i, 3), -q2(i, 2), q2(i, 1)] * j1 + ...
        [q1(i, 1), -q1(i, 2), q1(i, 3), q1(i, 4);
        q1(i, 2), q1(i, 1), -q1(i, 4), q1(i, 3);
        q1(i, 3), q1(i, 4), q1(i, 1), -q1(i, 2);
        q1(i, 4), -q1(i, 3), q1(i, 2), q1(i, 1)] * j2;
    jacob = [q3(i, 1), -q3(i, 2), -q3(i, 3), -q3(i, 4);
        q3(i, 2), q3(i, 1), q3(i, 4), -q3(i, 3);
        q3(i, 3), -q3(i, 4), q3(i, 1), q3(i, 2);
        q3(i, 4), q3(i, 3), -q3(i, 2), q3(i, 1)] * jm + ...
        [qm(i, 1), -qm(i, 2), qm(i, 3), qm(i, 4);
        qm(i, 2), qm(i, 1), -qm(i, 4), qm(i, 3);
        qm(i, 3), qm(i, 4), qm(i, 1), -qm(i, 2);
        qm(i, 4), -qm(i, 3), qm(i, 2), qm(i, 1)] * j3;
    grad_q(:, :, i) = grad_norm(:, :, i) * jacob;
end
end