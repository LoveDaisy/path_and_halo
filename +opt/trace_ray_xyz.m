function [ray_out, g_out] = trace_ray_xyz(ray_in, crystal, trace)
% INPUT
%    ray_in:        n*3, [longitude, latitude] in degree
%    crystal:       struct
%       .face_norm  m*3, [nx, ny, nz] face normals
%    trace:         struct
%       .n          k-vector, refract index after each face
%       .fid        k-vector, face id

face_norm = crystal.face_norm(trace.fid, :);
refract_n = opt.generate_trace_n(crystal, trace);
face_n = size(face_norm, 1);
ray_n = size(ray_in, 1);

% normalize
[ray_in, g_norm] = geo.normalize_vector(ray_in);

valid_idx = ray_in * face_norm(1, :)' < 0;

% refract on each face
curr_ray = ray_in;
g_out = g_norm;
for i = 1:face_n
    if i > 1
        valid_idx = valid_idx & (curr_ray * face_norm(i, :)' > 0);
    end
    n0 = refract_n(i);
    n1 = refract_n(i + 1);
    if n0 * n1 > 0
        [curr_ray, g_curr] = opt.refract(curr_ray, face_norm(i, :), abs(n0), abs(n1));
    else
        [curr_ray, g_curr] = opt.reflect(curr_ray, face_norm(i, :));
    end
    
    valid_idx = valid_idx & (~isnan(curr_ray(:, 1)));

    for j = 1:ray_n
        g_out(:, :, j) = g_curr(:, :, j) * g_out(:, :, j);
    end
end
ray_out = curr_ray;
ray_out(~valid_idx, :) = nan;
end

