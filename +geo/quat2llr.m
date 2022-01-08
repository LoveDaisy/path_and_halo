function llr = quat2llr(quat)
% Convert a rotation from quaternion to LLR representation.
% LLR is for [longitude, latitude, roll] and the corresponding rotation is: (for column vector, 
% rotation matrix is left-multiplied):
%   rotz(90 + longitude) * rotx(90 - latitude) * rotz(roll)
% See llr2mat() for detail.
%
% INPUT
%   quat:       n*4, quaternion
%
% OUTPUT
%   llr:        n*3, [longitude, latitude, roll] in degree

qz = quatrotate(quat, [0, 0, 1]);  % [cos(lat)cos(lon), cos(lat)sin(lon), sin(lat)]
latitude = asind(qz(:, 3));
longitude = atan2d(qz(:, 2), qz(:, 1));

qx = quatrotate(quat, [1, 0, 0]);  % [..., ..., cos(lat)sin(roll)]
qy = quatrotate(quat, [0, 1, 0]);  % [..., ..., cos(lat)cos(roll)]
roll = atan2d(qx(:, 3), qy(:, 3));
llr = [longitude, latitude, roll];
end