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

[r0, r0_ll, dr] = generate_healpix_grids(5);

r1 = nan(size(r0));
r2 = nan(size(r0));

entry_face_idx = 3;
exit_face_idx = 5;

%%
target_bending = 43;

res_store = find_bending_angle_contour(target_bending, ...
    face_norm([entry_face_idx, exit_face_idx], :), [1.31, 1]);

%%
figure(1); clf;
hold on;
for i = 1:length(res_store)
    plot(res_store{i}(:, 1), res_store{i}(:, 2), '-o');
    axis equal; axis tight;
    set(gca, 'xlim', [-180, 180], 'ylim', [-90, 90]);
end
box on;
