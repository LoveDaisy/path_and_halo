function [area_factor, vtx] = face_intersection_factor(vtx_entry, vtx_exit, ray_xyz)
% INPUT
%   vtx_entry:  entry face vertex
%   vtx_exit:   exit face vertex. it must be *CONVEX*
%   ray_xyz:    projection ray

e1 = vtx_exit(2, :) - vtx_exit(1, :);
e2 = vtx_exit(4, :) - vtx_exit(1, :);
e3 = cross(e1, e2);

uv0 = bsxfun(@minus, vtx_exit, vtx_exit(1, :)) / [e1; e2; e3];
uv0 = uv0(:, 1:2);
uvt = bsxfun(@minus, vtx_entry, vtx_exit(1, :)) / [e1; e2; -ray_xyz];
uvq = uvt(:, 1:2);

mask_vtx_num = size(uv0, 1);
q_vtx_num = size(uvq, 1);
uv = nan(q_vtx_num * 2, 2);

z_eps = 1e-8;
for i1 = 1:mask_vtx_num
    i2 = mod(i1, mask_vtx_num) + 1;
    i3 = mod(i2, mask_vtx_num) + 1;

    z0 = (uv0(i2, 1) - uv0(i1, 1)) * (uv0(i3, 2) - uv0(i2, 2)) - ...
        (uv0(i2, 2) - uv0(i1, 2)) * (uv0(i3, 1) - uv0(i2, 1));
    z2 = nan;
    uv_i = 1;
    for j1 = 1:q_vtx_num
        j2 = mod(j1, q_vtx_num) + 1;
        if ~isnan(z2)
            z1 = z2;
        else
            z1 = (uv0(i2, 1) - uv0(i1, 1)) * (uvq(j1, 2) - uv0(i1, 2)) - ...
                (uv0(i2, 2) - uv0(i1, 2)) * (uvq(j1, 1) - uv0(i1, 1));
        end
        z2 = (uv0(i2, 1) - uv0(i1, 1)) * (uvq(j2, 2) - uv0(i1, 2)) - ...
            (uv0(i2, 2) - uv0(i1, 2)) * (uvq(j2, 1) - uv0(i1, 1));

        if z0 * z1 > z_eps && z0 * z2 < z_eps
            % From inner to outer: record tow points
            uv(uv_i, :) = uvq(j1, :);
            ab =- (uv0(i1, :) - uvq(j1, :)) / [uv0(i2, :) - uv0(i1, :); uvq(j2, :) - uvq(j1, :)];
            uv(uv_i + 1, :) = uv0(i1, :) + ab(1) * (uv0(i2, :) - uv0(i1, :));
            uv_i = uv_i + 2;
        elseif z0 * z1 < -z_eps && z0 * z2 > -z_eps
            % From outer to inner: record one point
            ab =- (uv0(i1, :) - uvq(j1, :)) / [uv0(i2, :) - uv0(i1, :); uvq(j2, :) - uvq(j1, :)];
            uv(uv_i, :) = uv0(i1, :) + ab(1) * (uv0(i2, :) - uv0(i1, :));
            uv_i = uv_i + 1;
        elseif z0 * z1 > z_eps && z0 * z2 > z_eps
            % Two inner: record one point
            uv(uv_i, :) = uvq(j1, :);
            uv_i = uv_i + 1;
        elseif abs(z1) < z_eps
            uv(uv_i, :) = uvq(j1, :);
            uv_i = uv_i + 1;
        end
    end
    q_vtx_num = uv_i - 1;
    uvq = uv(1:q_vtx_num, :);
end
uv = uv(1:q_vtx_num, :);
vtx = uv * [vtx_exit(2, :) - vtx_exit(1, :); vtx_exit(4, :) - vtx_exit(1, :)];
vtx = bsxfun(@plus, vtx, vtx_exit(1, :));

if size(uv, 1) > 2
    area1 = sum(uv(:, 1) .* [uv(end, 2); uv(1:end - 1, 2)]) - ...
        sum([uv(end, 1); uv(1:end - 1, 1)] .* uv(:, 2));
    area0 = sum(uvt(:, 1) .* [uvt(end, 2); uvt(1:end - 1, 2)]) - ...
        sum([uvt(end, 1); uvt(1:end - 1, 1)] .* uvt(:, 2));
    area_factor = area1 / area0;
else
    area_factor = 0;
end
end
