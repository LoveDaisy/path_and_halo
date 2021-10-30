function q = find_rotation(origin, target, varargin)
% INPUT
%   origin:         n*6, [dx_in, dy_in, dz_in, dx_out, dy_out, dz_out]
%   target:         n*6
% PARAMETER (optional)
%   'eps':          scalar, default 1e-3, in degree
% OUTPUT
%   q:              n*4, quaternion

p = inputParser;
p.addParameter('eps', 1e-3);
p.parse(varargin{:});

num = size(origin, 1);

if p.Results.eps < inf
    bending_angle_origin = sum(origin(:, 1:3) .* origin(:, 4:6), 2) ./ ...
        sqrt(sum(origin(:, 1:3).^2, 2)) ./ sqrt(sum(origin(:, 4:6).^2, 2));
    bending_angle_target = sum(target(:, 1:3) .* target(:, 4:6), 2) ./ ...
        sqrt(sum(target(:, 1:3).^2, 2)) ./ sqrt(sum(target(:, 4:6).^2, 2));
    valid_idx = abs(bending_angle_origin - bending_angle_target) < p.Results.eps;
else
    valid_idx = true(num, 1);
end

q = zeros(num, 4);
if ~any(valid_idx)
    return;
end

origin = origin(valid_idx, :);
target = target(valid_idx, :);
num = size(origin, 1);

% find rotation between r0 and ray_in
origin_in = normalize_vector(origin(:, 1:3));
target_in = normalize_vector(target(:, 1:3));
q1_axis = zeros(num, 3);
for j = 1:num
    q1_axis(j, :) = cross(origin_in(j, :), target_in(j, :));
end
tmp_cos = sum(origin_in .* target_in, 2);
q1_theta = atan2d(sqrt(sum(q1_axis.^2, 2)), tmp_cos);
q1_axis = normalize_vector(q1_axis);
q1 = [cosd(-q1_theta/2), bsxfun(@times, sind(-q1_theta/2), q1_axis)];

% find rotation between rotated r2 and ray_out
origin_out = normalize_vector(origin(:, 4:6));
target_out = normalize_vector(target(:, 4:6));

origin_out = quatrotate(q1, origin_out);
origin_out = origin_out - bsxfun(@times, sum(origin_out .* target_in, 2), target_in);
origin_out = normalize_vector(origin_out);

target_out = target_out - bsxfun(@times, sum(target_out .* target_in, 2), target_in);
target_out = normalize_vector(target_out);
q2_axis = target_in;
q2_theta = acosd(sum(origin_out .* target_out, 2));
q2 = [cosd(-q2_theta/2), bsxfun(@times, sind(-q2_theta/2), q2_axis)];
q3 = [cosd(q2_theta/2), bsxfun(@times, sind(q2_theta/2), q2_axis)];
curr_idx = sum(quatrotate(q2, origin_out) .* target_out, 2) > ...
    sum(quatrotate(q3, origin_out) .* target_out, 2);
q2(~curr_idx, :) = q3(~curr_idx, :);

% total quaternion
q_all = quatmultiply(q1, q2);
q(valid_idx, :) = q_all;
end