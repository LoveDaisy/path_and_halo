function [x_contour, g_angle, y_val, jacobian] = ...
    find_bending_angle_contour(target_angle, face_norm, n, varargin)
% INPUT
%   target_angle:       scalar
%   face_norm:          n*3
%   n:                  n-vector

p = inputParser;
p.addParameter('eps', 1e-8);
p.addParameter('MaxIter', 100);
p.addParameter('GridLevel', 5);
p.addParameter('config', []);
p.parse(varargin{:});

if isempty(p.Results.config)
    config = generate_init_config(face_norm, n, p.Results.GridLevel);
else
    config = p.Results.config;
end
bending_angle_diff = abs(target_angle - config.bending_angle);
x_contour = {};
g_angle = {};
y_val = {};
jacobian = {};
if target_angle < config.bending_angle_min || target_angle > config.bending_angle_max
    return;
end


r_lim = config.dr * 2;

seeds_idx = bending_angle_diff < r_lim;
checked_idx = false(size(seeds_idx));

[~, min_idx] = min(bending_angle_diff);
start_p = config.r0_ll(min_idx, :);
i = 1;
while sum(seeds_idx & ~checked_idx) > 0
    [x_store_fwd, g_a_store_fwd, y_store_fwd, g_y_store_fwd, closed] = search_one_direction(start_p, ...
        target_angle, face_norm, n, 1, p.Results.MaxIter);
    valid_idx_fwd = ~isnan(x_store_fwd(:, 1));
    if ~closed
        [x_store_bck, g_a_store_bck, y_store_bck, g_y_store_bck, ~] = search_one_direction(start_p, ...
            target_angle, face_norm, n, -1, p.Results.MaxIter);
        valid_idx_bck = find(~isnan(x_store_bck(:, 1)));
        
        % filter out those identical to the other line
        d = distance_to_poly_line(x_store_bck(valid_idx_bck, 1:2), x_store_fwd(valid_idx_fwd, 1:2));
        valid_idx_bck = valid_idx_bck(d > r_lim * 2);
        
        valid_idx_bck = wrev(valid_idx_bck);
        curr_x = [x_store_bck(valid_idx_bck, 1:2); x_store_fwd(valid_idx_fwd, 1:2)];
        curr_g_a = [g_a_store_bck(valid_idx_bck, 1:2); g_a_store_fwd(valid_idx_fwd, 1:2)];
        curr_y = [y_store_bck(valid_idx_bck, 1:2); y_store_fwd(valid_idx_fwd, 1:2)];
        curr_g_y = cat(3, g_y_store_bck(:, :, valid_idx_bck), g_y_store_fwd(:, :, valid_idx_fwd));
    else
        curr_x = x_store_fwd(valid_idx_fwd, 1:2);
        curr_x = [curr_x; curr_x(1, :)];
        curr_g_a = g_a_store_fwd(valid_idx_fwd, 1:2);
        curr_g_a = [curr_g_a; curr_g_a(1, :)];
        curr_y = y_store_fwd(valid_idx_fwd, 1:2);
        curr_y = [curr_y; curr_y(1, :)];
        curr_g_y = g_y_store_fwd(:, :, valid_idx_fwd);
        curr_g_y = cat(3, curr_g_y, curr_g_y(:, :, 1));
    end
    
    if isempty(curr_x) || size(curr_x, 1) < 2
        break;
    end
    x_contour{i} = curr_x;
    g_angle{i} = curr_g_a;
    y_val{i} = curr_y;
    jacobian{i} = curr_g_y;
    i = i + 1;
    
    if ~closed
        % first filter out those very close to contour line
        next_idx = find(seeds_idx & ~checked_idx);
        d = distance_to_poly_line(config.r0_ll(next_idx, :), curr_x);
        checked_idx(next_idx(d < r_lim)) = true;
        next_idx = find(seeds_idx & ~checked_idx);
        
        % then find exact solution from other points
        for j = 1:length(next_idx)
            [x, a] = find_bending_angle_solution(config.r0_ll(next_idx(j), :), target_angle, ...
                face_norm, n, 'eps', r_lim * 0.1);
            d = distance_to_poly_line(x, curr_x);
            if abs(a - target_angle) > r_lim * 0.1 || d < r_lim
                checked_idx(next_idx(j)) = true;
            end
        end
        next_idx = find(seeds_idx & ~checked_idx);
        if isempty(next_idx)
            continue;
        end
        
        [~, min_idx] = min(bending_angle_diff(next_idx));
        start_p = config.r0_ll(next_idx(min_idx), :);
    else
        checked_idx(:) = true;
    end
end

dup_idx = false(length(x_contour), 1);
for i = 2:length(x_contour)
    curr_x = x_contour{i};
    curr_g_a = g_angle{i};
    curr_y = y_val{i};
    curr_g_y = jacobian{i};
    for j = 1:i-1
        prev_contour = x_contour{j};
        d = distance_to_poly_line(curr_x, prev_contour);
        valid_idx = d > r_lim * 2;
        
        curr_x = curr_x(valid_idx, :);
        curr_g_a = curr_g_a(valid_idx, :);
        curr_y = curr_y(valid_idx, :);
        curr_g_y = curr_g_y(:, :, valid_idx);
    end
    if isempty(curr_x) ||  size(curr_x, 1) < 2
        dup_idx(i) = true;
    end
    x_contour{i} = curr_x;
    g_angle{i} = curr_g_a;
    y_val{i} = curr_y;
    jacobian{i} = curr_g_y;
end
x_contour = x_contour(~dup_idx);
g_angle = g_angle(~dup_idx);
y_val = y_val(~dup_idx);
jacobian = jacobian(~dup_idx);
end


function [x_store, g_a_store, y_store, g_y_store, closed] = ...
    search_one_direction(x0, target_angle, face_norm, n, direction, res_num)
h = 1;
max_h = 5;
min_h = 0.1;
hh = -0.00;

[x, a, g_a, y, g_y] = find_bending_angle_solution(x0, target_angle, face_norm, n);

x_store = nan(res_num, 5);  % [x1, x2, angle]
x_store(1, :) = [x, a, h, 0];

y_store = nan(res_num, 2);
y_store(1, :) = y;

g_y_store = nan(2, 2, res_num);
g_y_store(:, :, 1) = g_y;

g_a_store = nan(res_num, 4);  % [original_gradient, normalized_co_gradient]
g_a_store(1, 1:2) = g_a;
g_a = g_a * [0, -1; 1, 0] * direction;
g_a_store(1, 3:4) = g_a / norm(g_a);

i = 2;
non_shrink_cnt = 0;
closed = false;

while i <= res_num
    if i <= 3
        x0 = x_store(i-1, 1:2) + h * g_a_store(i-1, 3:4);
    else
        % use 2nd order prediction
        tmp_t = cumsum(x_store(i-3:i-1, 5) * direction);
        tmp_x = x_store(i-3:i-1, 1:2);
        new_t = tmp_t(3) + h * direction;
        new_x = (new_t - tmp_t(1)) * (new_t - tmp_t(2)) * (tmp_t(1) - tmp_t(2)) * tmp_x(3, :) + ...
            (new_t - tmp_t(2)) * (new_t - tmp_t(3)) * (tmp_t(2) - tmp_t(3)) * tmp_x(1, :) + ...
            (new_t - tmp_t(3)) * (new_t - tmp_t(1)) * (tmp_t(3) - tmp_t(1)) * tmp_x(2, :);
        x0 = -new_x / ((tmp_t(1) - tmp_t(2)) * (tmp_t(2) - tmp_t(3)) * (tmp_t(3) - tmp_t(1)));
    end
    x0 = x0 + hh * h * g_a_store(i-1, 1:2) / norm(g_a_store(i-1, 1:2));
    [x, a, g_a, y, g_y] = find_bending_angle_solution(x0, target_angle, face_norm, n);
    if i >= 3 && dot(x - x_store(i-1, 1:2), g_a_store(i-1, 3:4)) / ...
            norm(x - x_store(i-1, 1:2)) < 0
        break;
    end
    
    if i > 3
        d = distance_to_poly_line(x, x_store(1:i-1, 1:2));
        if d < norm(x - x_store(i-1, 1:2)) * 0.5
            closed = true;
            break;
        end
    end
    
    g_a_store(i, 1:2) = g_a;
    g_a = g_a * [0, -1; 1, 0] * direction;
    g_a = g_a / norm(g_a);
    
    % a simple adaptive schedule
    if isnan(a)
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
    s = norm(x - x_store(i-1, 1:2));
    x_store(i, :) = [x, a, h, s];
    g_a_store(i, 3:4) = g_a;
    y_store(i, :) = y;
    g_y_store(:, :, i) = g_y;
    
    i = i + 1;
end
end