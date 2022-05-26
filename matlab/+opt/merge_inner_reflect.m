function [ray_out, jac] = merge_inner_reflect(ray_in, crystal, fid)
% Merge multiple inner reflect and reduce computation
%
% INPUT
%   ray_in:         n*3, [x, y, z], input ray directions. They may NOT be unit vectors.
%   crystal:        struct
%   fid:            k-vector, reflect faces
%
% OUTPUT
%   ray_out:        n*3
%   jac:            3*3*n, Jacobian, which takes ray_in as input and ray_out as output.

f_cnt = length(fid);
num = size(ray_in, 1);
need_jacobian = nargout == 2;

m = eye(3);
for i = 1:f_cnt
    m = (eye(3) - 2 * crystal.face_norm(fid(i), :)' * crystal.face_norm(fid(i), :)) * m;
end

ray_out = ray_in * m';
if need_jacobian
    jac = repmat(m, [1, 1, num]);
end
end
