function [mat, jac] = llr2mat(llr)
% Rotate matrix corresponding to llr-representation, i.e. [longitude, latitude, roll]
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
%   llr:        n*3, [longitude, latitude, roll] in degree
%
% OUTPUT
%   mat:        3*3*n, rotation matrix. It is for column vector, i.e. new_v = mat * v.
%   jac:        (3*3)*3*n, Jacobian, reshape first dimesion as matrix form

num = size(llr, 1);
c1 = cosd(llr(:, 1));
s1 = sind(llr(:, 1));
c2 = cosd(llr(:, 2));
s2 = sind(llr(:, 2));
c3 = cosd(llr(:, 3));
s3 = sind(llr(:, 3));

% The rotation matrix is:
%   -s1 * c3 - c1 * s2 * s3,  s1 * s3 - c1 * s2 * c3, c1 * c2
%    c1 * c3 - s1 * s2 * s3, -c1 * s3 - s1 * s2 * c3, s1 * c2
%    c2 * s3,                 c2 * c3,                s2

mat = zeros(3, 3, num);
for i = 1:num
    mat(:, :, i) = rotz(90 + llr(i, 1)) * rotx(90 - llr(i, 2)) * rotz(llr(i, 3));
end

j1 = [-c1 .* c3 + s1 .* s2 .* s3, c1 .* s3 + s1 .* s2 .* c3, -s1 .* c2, ...
        - s1 .* c3 - c1 .* s2 .* s3, s1 .* s3 - c1 .* s2 .* c3, c1 .* c2, ...
        zeros(num, 3)] * pi / 180;
j2 = [-c1 .* c2 .* s3, -c1 .* c2 .* c3, -c1 .* s2, ...
        - s1 .* c2 .* s3, -s1 .* c2 .* c3, -s1 .* s2, ...
        - s2 .* s3, -s2 .* c3, c2] * pi / 180;
j3 = [s1 .* s3 - c1 .* s2 .* c3, s1 .* c3 + c1 .* s2 .* s3, zeros(num, 1), ...
        - c1 .* s3 - s1 .* s2 .* c3, -c1 .* c3 + s1 .* s2 .* s3, zeros(num, 1), ...
        c2 .* c3, -c2 .* s3, zeros(num, 1)] * pi / 180;
jac = zeros(3, 3, 3, num);
for i = 1:num
    jac(:, :, 1, i) = reshape(j1(i, :), [3, 3])';
    jac(:, :, 2, i) = reshape(j2(i, :), [3, 3])';
    jac(:, :, 3, i) = reshape(j3(i, :), [3, 3])';
end
end
