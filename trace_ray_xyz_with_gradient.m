function [ray_out, bending_angle, g_out, g_angle] = ...
    trace_ray_xyz_with_gradient(ray_in, crystal, trace)
% INPUT
%    ray_in:        n*3, [longitude, latitude] in degree
%    crystal:       struct
%       .face_norm  m*3, [nx, ny, nz] face normals
%    trace:         struct
%       .n          k-vector, refract index after each face
%       .fid        k-vector, face id

face_norm = crystal.face_norm(trace.fid, :);
refract_n = [1; trace.n(:)];
face_n = size(face_norm, 1);
ray_n = size(ray_in, 1);

% normalize
[ray_in, g_norm] = normalize_vector_with_gradient(ray_in);

valid_idx = ray_in * face_norm(1, :)' < 0;

% refract on each face
curr_ray = ray_in;
g_out = g_norm;
for i = 1:face_n
    if i > 1
        valid_idx = valid_idx & (curr_ray * face_norm(i, :)' > 0);
    end
    [curr_ray, g_curr] = refract_with_gradient(curr_ray, face_norm(i, :), refract_n(i), refract_n(i+1));
    
    valid_idx = valid_idx & (~isnan(curr_ray(:, 1)));

    for j = 1:ray_n
        g_out(:, :, j) = g_curr(:, :, j) * g_out(:, :, j);
    end
end
ray_out = curr_ray;
ray_out(~valid_idx, :) = nan;

% find out bending angle
cos_angle = sum(curr_ray .* ray_in, 2);
g_cos = zeros(ray_n, 3);
for i = 1:ray_n
    g_cos(i, :) = ray_in(i, :) * g_out(:, :, i) + curr_ray(i, :) * g_norm(:, :, i);
end
bending_angle = acosd(cos_angle);
g_angle = bsxfun(@times, -1 ./ sqrt(1 - cos_angle.^2), g_cos) * 180 / pi;
end

