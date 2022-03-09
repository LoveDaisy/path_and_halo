clear; close all; clc;

load raypaths.mat

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

% SPC = sub parhelia circle
SPC_idx = [1, 2, 3, 4, 5, 6, 53, 54, 55, 56];
SPC_feature = rp_feature(SPC_idx, :);
SPC_reflect_mat = reflect_mat(:, :, SPC_idx);
SPC_rp_eigv_mat = rp_eigv_mat(:, :, SPC_idx);

% PC = parhelia circle
PC_idx = [10, 11, 12, 13, 14, 16, 17, 20, 22, 23, 24, 26, 27, 28, 29, 31, 33, 35, 36, ...
    37, 46, 58, 60, 61, 62, 63, 65, 66, 67, 68, 69, 70, 81, 84, 85];
PC_feature = rp_feature(PC_idx, :);
PC_reflect_mat = reflect_mat(:, :, PC_idx);
PC_rp_eigv_mat = rp_eigv_mat(:, :, PC_idx);

% SS = sub sun
SS_idx = [7, 8, 9, 15, 18, 19, 21, 25, 30, 32, 34];
SS_feature = rp_feature(SS_idx, :);
SS_reflect_mat = reflect_mat(:, :, SS_idx);
SS_rp_eigv_mat = rp_eigv_mat(:, :, SS_idx);

