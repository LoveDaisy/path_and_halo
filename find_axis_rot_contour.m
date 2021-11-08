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
    config = init_config_3d(face_norm, refract_n, sun_ll, p.Results.GridLevel);
else
    config = p.Results.config;
end

x_contour = {};
y_val = {};
jacobian = {};

target_xyz = ll2xyz_with_gradient(target_ll);
target_diff = acosd(config.out_xyz * target_xyz');

seeds_idx = target_diff < config.dr;
checked_idx = false(size(seeds_idx));
dup_dr = config.dr;
seeds_dr = config.dr;

[~, min_idx] = min(target_diff);
start_rot = config.axis_rot_store(min_idx, :);
checked_idx(min_idx) = true;

contour_i = 1;
while sum(seeds_idx & ~checked_idx) > 0
    [x_contour_fwd, y_val_fwd, jacobian_fwd, closed] = search_direction(start_rot, sun_ll, target_ll, ...
        face_norm, refract_n, 1, p.Results.MaxIter);
    valid_idx_fwd = find(~isnan(x_contour_fwd(:, 1)));
    if length(valid_idx_fwd) < 2
        valid_idx_fwd = [];
    end
    
    % check seeds
    for i = 2:length(valid_idx_fwd)
        tmp_idx = find(seeds_idx & ~checked_idx);
        [tmp_dist, tmp_vn, ~] = distance_to_poly_line(config.axis_rot_store(tmp_idx, :), ...
            x_contour_fwd(valid_idx_fwd(i-1:i), :));
        tmp_dy = tmp_vn * jacobian_fwd(:, :, i)';
        checked_idx(tmp_idx(tmp_dist < config.dr | ...
            (sum(tmp_dy.^2, 2) < config.dr^2))) = true;
    end

    if ~closed
        [x_contour_bck, y_val_bck, jacobian_bck, ~] = search_direction(start_rot, sun_ll, target_ll, ...
            face_norm, refract_n, -1, p.Results.MaxIter);
        valid_idx_bck = ~isnan(x_contour_bck(:, 1));
        valid_idx_bck = find(valid_idx_bck);
        if length(valid_idx_fwd) >= 2
            valid_idx_bck = valid_idx_bck(2:end);
        end

        % check seeds
        for i = 2:length(valid_idx_bck)
            tmp_idx = find(seeds_idx & ~checked_idx);
            [tmp_dist, tmp_vn, ~] = distance_to_poly_line(config.axis_rot_store(tmp_idx, :), ...
                x_contour_bck(valid_idx_bck(i-1:i), :));
            tmp_dy = tmp_vn * jacobian_bck(:, :, i)';
            checked_idx(tmp_idx(tmp_dist < config.dr | ...
                (sum(tmp_dy.^2, 2) < config.dr^2))) = true;
        end

        % filter out those identical to the other line
        num_before = length(valid_idx_bck);
        if ~isempty(valid_idx_bck) && length(valid_idx_fwd) > 1
            dist = distance_to_poly_line(x_contour_bck(valid_idx_bck, :), x_contour_fwd(valid_idx_fwd, :));
            valid_idx_bck = valid_idx_bck(dist > dup_dr);
        end
        if ~isempty(valid_idx_bck) && length(valid_idx_fwd) > 1
            tmp_x = x_contour_fwd(valid_idx_fwd, :);
            tmp_x(:, 3) = tmp_x(:, 3) + 360;
            dist = distance_to_poly_line(x_contour_bck(valid_idx_bck, :), tmp_x);
            valid_idx_bck = valid_idx_bck(dist > dup_dr);
        end
        if ~isempty(valid_idx_bck) && length(valid_idx_fwd) > 1
            tmp_x = x_contour_fwd(valid_idx_fwd, :);
            tmp_x(:, 3) = tmp_x(:, 3) - 360;
            dist = distance_to_poly_line(x_contour_bck(valid_idx_bck, :), tmp_x);
            valid_idx_bck = valid_idx_bck(dist > dup_dr);
        end
        num_after = length(valid_idx_bck);
        valid_idx_bck = wrev(valid_idx_bck);

        curr_x = [x_contour_bck(valid_idx_bck, :); x_contour_fwd(valid_idx_fwd, :)];
        curr_j = cat(3, jacobian_bck(:, :, valid_idx_bck), jacobian_fwd(:, :, valid_idx_fwd));
        curr_y = [y_val_bck(valid_idx_bck, :); y_val_fwd(valid_idx_fwd, :)];
        
        % entire contour is closed
        if num_before - num_after > 1 && norm(curr_x(1, :) - curr_x(end, :)) < dup_dr + ...
                max(norm(curr_x(2, :) - curr_x(1, :)), norm(curr_x(end, :) - curr_x(end-1, :)))
            curr_x = [curr_x; curr_x(1, :)];
            curr_j = cat(3, curr_j, curr_j(:, :, 1));
            curr_y = [curr_y; curr_y(1, :)];
        end
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
    
    % check seeds
    tmp_idx = find(seeds_idx & ~checked_idx);
    dist = distance_to_poly_line(config.axis_rot_store(tmp_idx, :), curr_x);
    checked_idx(tmp_idx(dist <= seeds_dr)) = true;
    tmp_idx = find(seeds_idx & ~checked_idx);
    for j = 1:length(tmp_idx)
        tmp_x0 = config.axis_rot_store(tmp_idx(j), :);
        [tmp_x, ~, ~] = find_solution(tmp_x0, sun_ll, target_ll, face_norm, refract_n);
        lon_offset = [-360, 1; 360, 1; 0, 1];
        lat_offset = [-180, -1; 180, -1; 0, 1];
        roll_offset = [-360, 1; 360, 1; 0, 1];
        comb_offset = zeros(27, 6);
        comb_i = 1;
        for lon_i = 1:3
            for lat_i = 1:3
                for roll_i = 1:3
                    comb_offset(comb_i, :) = ...
                        [lon_offset(lon_i, 1), lat_offset(lat_i, 1), roll_offset(roll_i, 1), ...
                        lon_offset(lon_i, 2), lat_offset(lat_i, 2), roll_offset(roll_i, 2)];
                    comb_i = comb_i + 1;
                end
            end
        end
        for comb_i = 1:27
            dist = distance_to_poly_line(tmp_x .* comb_offset(comb_i, 4:6) + ...
                comb_offset(comb_i, 1:3), curr_x);
            if dist <= seeds_dr
                checked_idx(tmp_idx(j)) = true;
                break;
            end
        end
    end
    tmp_idx = find(seeds_idx & ~checked_idx);
    
    [~, min_idx] = min(target_diff(tmp_idx));
    start_rot = config.axis_rot_store(tmp_idx(min_idx), :);
    
    % filter out identical to other lines up to a period
    dup = false;
    for j = 1:length(x_contour)
        tmp_x = curr_x;
        d1 = distance_to_poly_line(tmp_x, x_contour{j});
        tmp_x(:, 3) = curr_x(:, 3) + 360;
        d2 = distance_to_poly_line(tmp_x, x_contour{j});
        tmp_x(:, 3) = curr_x(:, 3) - 360;
        d3 = distance_to_poly_line(tmp_x, x_contour{j});
        if sum(min([d1, d2, d3], [], 2) < dup_dr) > length(d1) * 0.3
            dup = true;
            break;
        end
    end
    if dup
        continue;
    end
    
    x_contour{contour_i} = curr_x;
    y_val{contour_i} = curr_y;
    jacobian{contour_i} = curr_j;
    contour_i = contour_i + 1;
end
end


function [x_contour, y_val, jacobian, closed] = search_direction(rot0, sun_ll, target_ll, face_norm, refract_n, ...
    direction, num)
h = 1;
max_h = 5;
min_h = 0.3;

x_contour = nan(num, 3);
y_val = nan(num, 2);
jacobian = nan(2, 3, num);
co_g_store = nan(num, 3);
ds_store = nan(num, 1);
closed = false;

[rot, out_ll, j_rot] = find_solution(rot0, sun_ll, target_ll, face_norm, refract_n);
co_gradient = cross(j_rot(1, :), j_rot(2, :));
co_gradient = co_gradient / norm(co_gradient) * direction;

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
        tmp_t = cumsum(ds_store(i-3:i-1));
        tmp_x = x_contour(i-3:i-1, :);
        new_t = tmp_t(3) + h;
        new_x = (new_t - tmp_t(1)) * (new_t - tmp_t(2)) * (tmp_t(1) - tmp_t(2)) * tmp_x(3, :) + ...
            (new_t - tmp_t(2)) * (new_t - tmp_t(3)) * (tmp_t(2) - tmp_t(3)) * tmp_x(1, :) + ...
            (new_t - tmp_t(3)) * (new_t - tmp_t(1)) * (tmp_t(3) - tmp_t(1)) * tmp_x(2, :);
        x0 = -new_x / ((tmp_t(1) - tmp_t(2)) * (tmp_t(2) - tmp_t(3)) * (tmp_t(3) - tmp_t(1)));
    end
    [x, out_ll, j_rot] = find_solution(x0, sun_ll, target_ll, face_norm, refract_n);
    co_gradient = cross(j_rot(1, :), j_rot(2, :));
    co_gradient = co_gradient / norm(co_gradient) * direction;
    co_g_store(i, :) = co_gradient;

    if i >= 3 && dot(x - x_contour(i-1, :), co_g_store(i, :)) / ...
            norm(x - x_contour(i-1, :)) < 0
        break;
    end

    if i > 3
        [d, ~, idx_t] = distance_to_poly_line(x, x_contour(1:i-1, :));
        if d < norm(x - x_contour(i-1, :)) * 0.2
            closed = true;
            if idx_t(2) > 1 - 1e-6 || idx_t(2) < 1e-6
                x_contour(i, :) = x;
                jacobian(:, :, i) = j_rot;
                y_val(i, :) = out_ll;
            end
            break;
        end
    end

    % a simple adaptive schedule
    % 1. estimate curvature
    theta = acos(dot(x_contour(i-1, :), x) / norm(x_contour(i-1, :)) / norm(x));
    rho = sqrt(sum((x_contour(i-1, :) - x).^2, 2)) / theta;
    
    % 2. shrink or extend
    if any(isnan(x)) || abs(norm(x - x_contour(i-1, :)) - h) / h > 0.7 || theta > 2 * pi / 180
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
        h = max(min(h, rho * 0.02), min_h);
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
p.addParameter('MaxIter', 10);
p.parse(varargin{:});

max_step = 50;

out_ll = nan(1, 2);
j_rot = nan(2, 3, 1);
x = nan(1, 3);

if any(isnan(rot0)) || any(isinf(rot0))
    return;
end

[out_ll, j_rot] = crystal_system_with_gradient(rot0, sun_ll, face_norm, refract_n);

if any(isnan(out_ll))
    return;
end

dy = target_ll - out_ll;
iter_num = 1;
x = rot0;
while norm(dy) > p.Results.eps && iter_num < p.Results.MaxIter
    dx = dy / j_rot';
    dx = min(norm(dx), max_step) * dx / norm(dx);

    % then linear search
    out_ll = nan;
    alpha = 2;
    while (any(isnan(out_ll)) || norm(target_ll - out_ll) > norm(dy) * 0.8) && alpha > 0.1
        alpha = alpha / 2;
        [out_ll, j_rot] = crystal_system_with_gradient(x + dx * alpha, sun_ll, face_norm, refract_n);
    end

    x = x + dx * alpha;
    dy = target_ll - out_ll;
    iter_num = iter_num + 1;
end

if any(isnan(out_ll)) || norm(dy) > p.Results.eps
    out_ll = nan(1, 2);
    j_rot = nan(2, 3, 1);
    x = nan(1, 3);
end
end