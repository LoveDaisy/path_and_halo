clear; close all; clc;

face_norm = [0, 0, 1;     % face 1
    0, 0, -1;             % face 2
    -sqrt(3)/2, 1/2, 0;   % face 3
    0, 1, 0;              % face 4
    sqrt(3)/2, 1/2, 0;    % face 5
    sqrt(3)/2, -1/2, 0;   % face 6
    0, -1, 0;             % face 7
    -sqrt(3)/2, -1/2, 0]; % face 8

n = [1.31; 1];

ray_in = [20, -35; 20.1, -35; 20, -34.95; -10, -10; -10.1, -10];
[ray_out, bending_angle, g_ray, g_angle] = bending_angle_with_gradient(ray_in, face_norm([3, 5], :), n);
ray_out_xyz = [cosd(ray_out(:, 2)) .* cosd(ray_out(:, 1)), ...
    cosd(ray_out(:, 2)) .* sind(ray_out(:, 1)), ...
    sind(ray_out(:, 2))];

