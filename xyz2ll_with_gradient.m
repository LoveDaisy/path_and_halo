function [ll, g] = xyz2ll_with_gradient(xyz)
ll = [atan2d(xyz(:, 2), xyz(:, 1)), ...
    asind(xyz(:, 3) ./ sqrt(sum(xyz.^2, 2)))];

m = size(xyz, 1);
g1 = [-xyz(:, 2) ./ sum(xyz(:, 1:2).^2, 2), ...
    -xyz(:, 1) .* xyz(:, 3) ./ sqrt(sum(xyz(:, 1:2).^2, 2)) ./ sum(xyz.^2, 2)];  % size m*2
g2 = [xyz(:, 1) ./ sum(xyz(:, 1:2).^2, 2), ...
    -xyz(:, 2) .* xyz(:, 3) ./ sqrt(sum(xyz(:, 1:2).^2, 2)) ./ sum(xyz.^2, 2)];  % size m*2
g3 = [zeros(m, 1), ...
    sqrt(sum(xyz(:, 1:2).^2, 2)) ./ sum(xyz.^2, 2)];
g = cat(3, g1, g2, g3) * 180 / pi;  % size m*2*3
g = permute(g, [2, 3, 1]);
end