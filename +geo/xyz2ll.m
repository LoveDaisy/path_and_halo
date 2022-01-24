function [ll, jac] = xyz2ll(xyz)
% Convert Cartesian coordinate [x, y, z] into spherical coordinate [longitude, latidue]
% Input vector [x, y, z] may not be unit vector.
%
% INPUT
%   xyz:        n*3, [x, y, z]
%
% OUTPUT
%   ll:         n*2, [longitude, latitude], in degree
%   jac:        2*3*n, Jacobian for this transformation

need_jacobian = nargout == 2;

ll = [atan2d(xyz(:, 2), xyz(:, 1)), ...
        asind(xyz(:, 3) ./ sqrt(sum(xyz.^2, 2)))];

if need_jacobian
    m = size(xyz, 1);
    g1 = [-xyz(:, 2) ./ sum(xyz(:, 1:2).^2, 2), ...
            -xyz(:, 1) .* xyz(:, 3) ./ sqrt(sum(xyz(:, 1:2).^2, 2)) ./ sum(xyz.^2, 2)]; % size m*2
    g2 = [xyz(:, 1) ./ sum(xyz(:, 1:2).^2, 2), ...
            -xyz(:, 2) .* xyz(:, 3) ./ sqrt(sum(xyz(:, 1:2).^2, 2)) ./ sum(xyz.^2, 2)]; % size m*2
    g3 = [zeros(m, 1), ...
            sqrt(sum(xyz(:, 1:2).^2, 2)) ./ sum(xyz.^2, 2)];
    jac = cat(3, g1, g2, g3) * 180 / pi; % size m*2*3
    jac = permute(jac, [2, 3, 1]);
end
end
