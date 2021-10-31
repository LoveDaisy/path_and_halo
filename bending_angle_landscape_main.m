clear; close all; clc;

face_norm = [0, 0, 1;     % face 1
    0, 0, -1;             % face 2
    -sqrt(3)/2, 1/2, 0;   % face 3
    0, 1, 0;              % face 4
    sqrt(3)/2, 1/2, 0;    % face 5
    sqrt(3)/2, -1/2, 0;   % face 6
    0, -1, 0;             % face 7
    -sqrt(3)/2, -1/2, 0]; % face 8
n = 1.31;

side_level = 5;
[r0, r0_ll, dr] = generate_healpix_grids(side_level);

r1 = nan(size(r0));
r2 = nan(size(r0));

entry_face_idx = 3;
exit_face_idx = 5;
entry_norm = face_norm(entry_face_idx, :);
exit_norm = face_norm(exit_face_idx, :);

valid_idx = r0 * entry_norm' < 0;
r1(valid_idx, :) = refract_with_gradient(r0(valid_idx, :), entry_norm, 1, n);
valid_idx = valid_idx & (r1 * exit_norm' > 0);
r2(valid_idx, :) = refract_with_gradient(r1(valid_idx, :), exit_norm, n, 1);
valid_idx = valid_idx & sum(r2.^2, 2) > 1e-4;
r2(~valid_idx, :) = nan;

r2_ll = xyz2ll_with_gradient(r2);

bending_angle = acosd(sum(r0 .* r2, 2));
bending_angle_max = max(bending_angle);
bending_angle_min = min(bending_angle);

target_bending = 43;

figure(1); clf;
hold on;
scatter(r0_ll(valid_idx, 1), r0_ll(valid_idx, 2), 20, bending_angle(valid_idx), 'fill');
axis equal; axis tight;
box on;
set(gca, 'xlim', [-180, 180], 'ylim', [-90, 90]);

idx = abs(target_bending - bending_angle) < dr;
plot(r0_ll(idx, 1), r0_ll(idx, 2), 'rs');