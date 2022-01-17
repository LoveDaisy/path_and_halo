function [ray_out, jac_out] = trace_ray_direction(ray_in, crystal, trace)
% Trace ray through a crystal. It is NOT a complete ray tracing, but only count ray direction.
%
% INPUT
%    ray_in:        n*3, [x, y, z], they may NOT be unit vectors
%    crystal:       struct
%    trace:         struct
%
% OUTPUT
%   ray_out:        n*3, [x, y, z], they are unit vectors
%   jac_out:        3*3*n, Jacobian, input is ray_in and output is ray_out

face_norm = crystal.face_norm(trace.fid, :);
refract_n = opt.generate_trace_n(crystal, trace);
face_cnt = size(face_norm, 1);
ray_cnt = size(ray_in, 1);

valid_idx = ray_in * face_norm(1, :)' < 0;

jac_out = zeros(3, 3, ray_cnt);
for i = 1:ray_cnt
    jac_out(:, :, i) = eye(3);
end

% refract/reflect on each face
curr_ray = ray_in;
i = 1;
while i <= face_cnt
    if i > 1
        valid_idx = valid_idx & (curr_ray * face_norm(i, :)' > 0);
    end
    n0 = refract_n(i);
    n1 = refract_n(i + 1);

    if face_cnt > 2 && i > 1 && i < face_cnt
        [curr_ray, jac_curr] = opt.merge_inner_reflect(curr_ray, crystal, trace.fid(2:end - 1));
        i = face_cnt - 1;
    elseif n0 * n1 > 0
        [curr_ray, jac_curr] = opt.refract(curr_ray, face_norm(i, :), abs(n0), abs(n1));
    else
        [curr_ray, jac_curr] = opt.reflect(curr_ray, face_norm(i, :));
    end

    valid_idx = valid_idx & (~isnan(curr_ray(:, 1)));

    for j = 1:ray_cnt
        jac_out(:, :, j) = jac_curr(:, :, j) * jac_out(:, :, j);
    end
    i = i + 1;
end
ray_out = curr_ray;
ray_out(~valid_idx, :) = nan;
end
