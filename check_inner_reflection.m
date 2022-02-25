clear; close all; clc;

load all_raypaths_5.mat

num = size(possible_raypaths, 1);
rp_len = sum(possible_raypaths > 0, 2);

reflect_mat = nan(3, 3, num);  % output_ray = mat * input_ray;
normal_mat = nan(2, 3, num);   % [n1 * mat'; n2]
feature_mat = nan(num, 7);
ev_real = nan(num, 3);
ev_imag = nan(num, 3);
for i = 1:num
    if rp_len(i) <= 2
        continue;
    end
    
    t = eye(3);
    for j = 2:rp_len(i)-1
        n = crystal.face_norm(possible_raypaths(i, j), :);
        t = (eye(3) - 2 * (n' * n)) * t;
    end
    [t_V, t_D] = eig(t);
    t_d = diag(t_D);
    reflect_mat(:, :, i) = t;
    ev_real(i, :) = real(t_d');
    ev_imag(i, :) = imag(t_d');
    
    n1 = crystal.face_norm(possible_raypaths(i, 1), :);
    n2 = crystal.face_norm(possible_raypaths(i, rp_len(i)), :);
    normal_mat(:, :, i) = [n1 * t'; n2];

    feature_mat(i, :) = [dot(n1 * t', n2), ev_real(i, :), ev_imag(i, :)];
end
feature_mat = round(feature_mat, 6);
[~, ia, ic] = unique(feature_mat, 'rows');

[~, is] = sort(ic);
label_rp = [ic(is), possible_raypaths(is, :)];

idx = rp_len(is) == 5;
curr_label_rp = label_rp(idx, :);
curr_feature = feature_mat(is(idx), :);
curr_reflect_mat = reflect_mat(:, :, is(idx));
