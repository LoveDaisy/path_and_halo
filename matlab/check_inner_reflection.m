clear; close all; clc;

load column_raypaths_6.mat

num = size(possible_raypaths, 1);
rp_len = sum(possible_raypaths > 0, 2);

reflect_mat = nan(3, 3, num);  % output_ray = mat * input_ray;
rp_feature = nan(num, 7);
rp_eigv_mat = nan(3, 6, num);
for i = 1:num
    if rp_len(i) < 2
        n = crystal.face_norm(possible_raypaths(i, 1), :);
        t = (eye(3) - 2 * (n' * n));
    else
        t = eye(3);
        for j = 2:rp_len(i)-1
            n = crystal.face_norm(possible_raypaths(i, j), :);
            t = (eye(3) - 2 * (n' * n)) * t;
        end
    end
    
    [t_V, t_D] = eig(t);
    t_d = diag(t_D);
    reflect_mat(:, :, i) = t;
    
    n1 = crystal.face_norm(possible_raypaths(i, 1), :);
    n2 = crystal.face_norm(possible_raypaths(i, rp_len(i)), :);

    rp_feature(i, :) = [dot(n1 * t', n2), real(t_d'), imag(t_d')];
    rp_eigv_mat(:, 1:3, i) = real(t_V);
    rp_eigv_mat(:, 4:6, i) = imag(t_V);
end
rp_feature = round(rp_feature, 6);
[~, ia, ic] = unique(rp_feature, 'rows');

[~, is] = sort(ic);
label_rp = [ic(is), possible_raypaths(is, :)];

rp_feature = rp_feature(is, :);
reflect_mat = reflect_mat(:, :, is);
rp_eigv_mat = rp_eigv_mat(:, :, is);

