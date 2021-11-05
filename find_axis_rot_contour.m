function [x_contour, y_val, jacobian] = find_axis_rot_contour(sun_ll, target_ll, ...
    face_norm, refract_n, varargin)
% INPUT
%   sun_ll:             2-vector, [longitude, latitude] in degree
%   target_ll:          2-vector, [longitude, latitude] in degree
%   face_norm:          k*3
%   refract_n:          k-vector

p = inputParser;
p.addParameter('eps', 1e-8);
p.addParameter('MaxIter', 100);
p.addParameter('GridLevel', 5);
p.addParameter('config', []);
p.parse(varargin{:});

if isempty(p.Results.config)
    config = init_config_3d(face_norm, n, sun_ll, p.Results.GridLevel);
else
    config = p.Results.config;
end

target_xyz = ll2xyz_with_gradient(target_ll);
target_diff = acosd(config.out_xyz * target_xyz');

seeds_idx = target_diff < config.dr;
checked_idx = false(size(seeds_idx));

[~, min_idx] = min(target_diff);
start_rot = config.axis_rot_store(min_idx, :);
contour_i = 1;
while sum(seeds_idx & ~checked_idx) > 0
    [x_contour_fwd, y_val_fwd, jacobian_fwd, closed] = search_one_direction(start_rot, sun_ll, target_ll, ...
        face_norm, refract_n, 1, num);
    valid_idx_fwd = ~isnan(x_contour_fwd(:, 1));

    if ~closed
        [x_contour_bck, y_val_bck, jacobian_bck, ~] = search_one_direction(start_rot, sun_ll, target_ll, ...
            face_norm, refract_n, 1, num);
        valid_idx_bck = ~isnan(x_contour_bck(:, 1));

        % filter out those identical to the other line
        d = distance_to_poly_line(x_contour_bck(valid_idx_bck, :), x_contour_fwd(valid_idx_fwd, :));
        valid_idx_bck = valid_idx_bck(d > config.dr * 2);
        
        valid_idx_bck = wrev(valid_idx_bck);
        curr_x = [x_contour_bck(valid_idx_bck, :); x_contour_fwd(valid_idx_fwd, :)];
        curr_j = cat(3, jacobian_bck(:, :, valid_idx_bck), jacobian_fwd(:, :, valid_idx_fwd));
        curr_y = [y_val_bck(valid_idx_bck, :); y_val_fwd(valid_idx_fwd, :)];
    else
        curr_x = x_contour_fwd(valid_idx_fwd, :);
        curr_x = [curr_x; curr_x(1, :)];
        curr_j = jacobian_fwd(:, :, valid_idx_fwd);
        curr_j = cat(3, curr_j, curr_j(:, :, 1));
        curr_y = y_val_fwd(valid_idx_fwd, :);
        curr_y = [curr_y; curr_y(1, :)];
    end

    if isempty(curr_x) || size(curr_x, 1) < 2
        break;
    end
    x_contour{contour_i} = curr_x;
    y_val{contour_i} = curr_y;
    jacobian{i} = curr_j;
    contour_i = contour_i + 1;
end
end


function [x_contour, y_val, jacobian, closed] = search_one_direction(rot0, sun_ll, target_ll, face_norm, refract_n, ...
    direction, num)
h = 1;
max_h = 5;
min_h = 0.1;

x_contour = nan(num, 3);
y_val = nan(num, 2);
jacobian = nan(2, 3, num);
co_g_store = nan(num, 3);
ds_store = nan(num, 1);
closed = false;

[rot, out_ll, j_rot] = find_solution(rot0, sun_ll, target_ll, face_norm, refract_n);
co_gradient = cross(j_rot(1, :), j_rot(2, :));
co_gradient = co_gradient' / norm(co_gradient) * direction;

x_contour(1, :) = rot;
y_val(1, :) = out_ll;
jacobian(:, :, 1) = j_rot;
co_g_store(1, :) = co_gradient;
ds_store(1) = 0;

i = 2;
non_shrink_cnt = 0;
while i <= num
    if i <= 3
        x0 = x_contour(i-1, :) + h * co_gradient;
    else
        tmp_t = cumsum(ds_store(i-3:i-1) * direction);
        tmp_x = x_contour(i-3:i-1, :);
        new_t = tmp_t(3) + h * direction;
        new_x = (new_t - tmp_t(1)) * (new_t - tmp_t(2)) * (tmp_t(1) - tmp_t(2)) * tmp_x(3, :) + ...
            (new_t - tmp_t(2)) * (new_t - tmp_t(3)) * (tmp_t(2) - tmp_t(3)) * tmp_x(1, :) + ...
            (new_t - tmp_t(3)) * (new_t - tmp_t(1)) * (tmp_t(3) - tmp_t(1)) * tmp_x(2, :);
        x0 = -new_x / ((tmp_t(1) - tmp_t(2)) * (tmp_t(2) - tmp_t(3)) * (tmp_t(3) - tmp_t(1)));
    end
    [x, out_ll, j_rot] = find_solution(x0, sun_ll, target_ll, face_norm, n);
    co_gradient = cross(j_rot(1, :), j_rot(2, :));
    co_gradient = co_gradient' / norm(co_gradient) * direction;

    if i >= 3 && dot(x - x_contour(i-1, :), co_g_store(i, :)) / ...
            norm(x - x_contour(i-1, :)) < 0
        break;
    end

    if i > 3
        d = distance_to_poly_line(x, x_contour(1:i-1, :));
        if d < norm(x - x_contour(i-1, :)) * 0.5
            closed = true;
            break;
        end
    end

    % a simple adaptive schedule
    if any(isnan(x))
        non_shrink_cnt = 0;
        h = h / 2;
        if h > min_h
            continue;
        else
            break;
        end
    elseif h < max_h
        non_shrink_cnt = non_shrink_cnt + 1;
        if  non_shrink_cnt > 1
            h = h * 1.2;
        end
    end
    ds_store(i) = norm(x - x_contour(i-1, :));
    x_contour(i, :) = x;
    jacobian(:, :, i) = j_rot;
    y_val(i, :) = out_ll;

    i = i + 1;
end
end


function [x, out_ll, j_rot] = find_solution(rot0, sun_ll, target_ll, face_norm, refract_n, varargin)
% INPUT
%   rot0:           1*3
%   sun_ll:         1*2
%   target_ll:      1*2
%   face_norm:      k*3
%   refract_n:      k*1

p = inputParser;
p.addParameter('eps', 1e-8);
p.addParameter('MaxIter', 5);
p.parse(varargin{:});

max_step = 50;

num = size(rot0, 1);
[out_ll, j_rot] = crystal_system_with_gradient(rot0, sun_ll, face_norm, refract_n);

if any(isnan(out_ll))
    out_ll = nan(num, 2);
    j_rot = nan(2, 3, num);
    return;
end

dy = target_ll - out_ll;
iter_num = 1;
x = rot0;
while norm(dy) > p.Results.eps && iter_num < p.Results.MaxIter
    dx = j_rot' * (j_rot * j_rot') \ dy;
    dx = min(norm(dx), max_step) * dx / norm(dx);

    % then linear search
    out_ll = nan;
    alpha = 2;
    while (any(isnan(out_ll)) || norm(target_ll - out_ll) > norm(dy) * 0.8) && alpha > 0.1
        alpha = alpha / 2;
        [out_ll, ~] = crystal_system_with_gradient(x + dx * alpha, face_norm, n);
    end

    x = x + dx;
    dy = target_ll - out_ll;
    iter_num = iter_num + 1;
end

if norm(dy) > p.Results.eps
    out_ll = nan(num, 2);
    j_rot = nan(2, 3, num);
end
end