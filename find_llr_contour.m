function [x_contour, y_val, jacobian] = find_llr_contour(sun_ll, target_ll, crystal, trace, varargin)
% INPUT
%   sun_ll:             2-vector, [longitude, latitude] in degree
%   target_ll:          2-vector, [longitude, latitude] in degree
%   crystal:
%   trace:

p = inputParser;
p.addParameter('eps', 1e-8);
p.addParameter('MaxIter', 100);
p.addParameter('GridLevel', 5);
p.addParameter('config', []);
p.addParameter('OutputSpace', 'llr');
p.parse(varargin{:});

if isempty(p.Results.config)
    config = init_config_3d(crystal, trace, sun_ll, p.Results.GridLevel);
else
    config = p.Results.config;
end

x_contour = {};
y_val = {};
jacobian = {};

target_xyz = geo.ll2xyz(target_ll);
target_diff = acosd(config.out_xyz * target_xyz');

dup_dr = config.dr;
seeds_dr = config.dr * 2;
seeds_idx = target_diff < seeds_dr;
checked_idx = false(size(target_diff));

[min_diff, min_idx] = min(target_diff(seeds_idx));
tmp_idx = find(seeds_idx);
start_x = config.axis_llr_store(tmp_idx(min_idx), :);
checked_idx(tmp_idx(min_idx)) = true;

contour_i = 1;
while sum(seeds_idx & ~checked_idx) > 0
    [x_contour_fwd, y_val_fwd, jacobian_fwd, closed] = search_direction(start_x, sun_ll, target_ll, ...
        crystal, trace, 1, p.Results.MaxIter);
    valid_idx_fwd = find(~isnan(x_contour_fwd(:, 1)));
    if length(valid_idx_fwd) < 2
        valid_idx_fwd = [];
    end

    if ~closed
        [x_contour_bck, y_val_bck, jacobian_bck, ~] = search_direction(start_x, sun_ll, target_ll, ...
            crystal, trace, -1, p.Results.MaxIter);
        valid_idx_bck = ~isnan(x_contour_bck(:, 1));
        valid_idx_bck = find(valid_idx_bck);
        if length(valid_idx_fwd) >= 2
            valid_idx_bck = valid_idx_bck(2:end);
        end

        % filter out those identical to the other line
        num_before = length(valid_idx_bck);
        if ~isempty(valid_idx_bck) && length(valid_idx_fwd) > 1
            [dist, ~, idx_t] = distance_to_poly_line(x_contour_bck(valid_idx_bck, :), ...
                x_contour_fwd(valid_idx_fwd, :));
            valid_idx_bck = valid_idx_bck(dist > dup_dr | ...
                (idx_t(:, 2) <= 1e-8 | idx_t(:, 2) >= 1-1e-8));
        end
        if ~isempty(valid_idx_bck) && length(valid_idx_fwd) > 1
            tmp_x = x_contour_fwd(valid_idx_fwd, :);
            tmp_x(:, 3) = tmp_x(:, 3) + 360;
            [dist, ~, idx_t] = distance_to_poly_line(x_contour_bck(valid_idx_bck, :), tmp_x);
            valid_idx_bck = valid_idx_bck(dist > dup_dr | ...
                (idx_t(:, 2) <= 1e-8 | idx_t(:, 2) >= 1-1e-8));
        end
        if ~isempty(valid_idx_bck) && length(valid_idx_fwd) > 1
            tmp_x = x_contour_fwd(valid_idx_fwd, :);
            tmp_x(:, 3) = tmp_x(:, 3) - 360;
            [dist, ~, idx_t] = distance_to_poly_line(x_contour_bck(valid_idx_bck, :), tmp_x);
            valid_idx_bck = valid_idx_bck(dist > dup_dr | ...
                (idx_t(:, 2) <= 1e-8 | idx_t(:, 2) >= 1-1e-8));
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
        checked_idx(target_diff >= min_diff & seeds_idx) = true;
        tmp_idx = find(seeds_idx & ~checked_idx);
        [min_diff, min_idx] = min(target_diff(tmp_idx));
        start_x = config.axis_llr_store(tmp_idx(min_idx), :);
        continue;
    end
    
    % check seeds
    tmp_idx = find(seeds_idx & ~checked_idx);
    [tmp_idist, tmp_odist] = input_output_distance(config.axis_llr_store(tmp_idx, :), curr_x, curr_j);
    tmp_dist_checked = tmp_odist <= seeds_dr | tmp_idist <= seeds_dr;
    checked_idx(tmp_idx) = tmp_dist_checked;
    tmp_idx = find(seeds_idx & ~checked_idx);
    for j = 1:length(tmp_idx)
        if checked_idx(tmp_idx(j))
            continue;
        end
        tmp_x0 = config.axis_llr_store(tmp_idx(j), :);
        [tmp_x, ~, ~] = find_solution(tmp_x0, sun_ll, target_ll, crystal, trace, ...
            'eps', config.dr * 0.25);
        if any(isnan(tmp_x)) && target_diff(tmp_idx(j)) > seeds_dr
            checked_idx(target_diff > target_diff(tmp_idx(j)) & seeds_idx) = true;
            continue;
        end

        [tmp_idist1, tmp_odist1] = input_output_distance(tmp_x, curr_x, curr_j);
        if tmp_idist1 <= seeds_dr || tmp_odist1 <= seeds_dr
            checked_idx(tmp_idx(j)) = true;
        end
    end
    tmp_idx = find(seeds_idx & ~checked_idx);
    
    [min_diff, min_idx] = min(target_diff(tmp_idx));
    start_x = config.axis_llr_store(tmp_idx(min_idx), :);
    checked_idx(tmp_idx(min_idx)) = true;
    
    % filter out identical to other lines up to a period
    dup_idx = false(size(curr_x, 1), 1);
    for j = 1:length(x_contour)
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
            tmp_x = bsxfun(@plus, bsxfun(@times, curr_x, comb_offset(comb_i, 4:6)), comb_offset(comb_i, 1:3));
            tmp_d = distance_to_poly_line(tmp_x, x_contour{j});
            dup_idx = dup_idx | tmp_d < dup_dr;
        end
    end
    if sum(dup_idx) / size(curr_x, 1) > 0.3
        continue;
    end
    
    x_contour{contour_i} = curr_x;
    y_val{contour_i} = curr_y;
    jacobian{contour_i} = curr_j;
    contour_i = contour_i + 1;
end
end


function [x_contour, y_val, jacobian, closed] = search_direction(rot0, sun_ll, target_ll, crystal, trace, ...
    direction, num)
h = 3;
max_h = 15;
min_h = 0.3;

x_contour = nan(num, 3);
y_val = nan(num, 2);
jacobian = nan(2, 3, num);
co_g_store = nan(num, 3);
ds_store = nan(num, 1);
closed = false;

[rot, out_ll, j_rot] = find_solution(rot0, sun_ll, target_ll, crystal, trace);
if any(isnan(out_ll))
    return;
end
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
    [x, out_ll, j_rot] = find_solution(x0, sun_ll, target_ll, crystal, trace);
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
    if i > 2
        dx1 = x - x_contour(i-1, :);
        dx0 = x_contour(i-1, :) - x_contour(i-2, :);
        theta = acos(dot(dx1, dx0) / norm(dx1) / norm(dx0));
        rho = norm(dx0) / theta;
    else
        dx1 = h;
        theta = 0;
        rho = inf;
    end
    
    % 2. shrink or extend
    if any(isnan(x)) || abs(norm(dx1) - h) / h > 0.7 || theta > 40 * pi / 180
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
        h = max(min(h, rho * .4), min_h);
    end
    ds_store(i) = norm(x - x_contour(i-1, :));
    x_contour(i, :) = x;
    jacobian(:, :, i) = j_rot;
    y_val(i, :) = out_ll;

    i = i + 1;
end
end


function [x, out_ll, j_rot] = find_solution(rot0, sun_ll, target_ll, crystal, trace, varargin)
% INPUT
%   rot0:           1*3
%   sun_ll:         1*2
%   target_ll:      1*2
%   crystal:
%   trace:

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

[out_ll, j_rot] = opt.crystal_system(rot0, sun_ll, crystal, trace);

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
        [out_ll, j_rot] = opt.crystal_system(x + dx * alpha, sun_ll, crystal, trace);
    end
    
    if norm(target_ll - out_ll) > norm(dy)
        dy = target_ll - out_ll;
        break;
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


function [dist, vn, idx_t] = distance_to_poly_line(qxy, xy)
% INPUT
%   qxy:    n*2
%   xy:     m*2

n = size(qxy, 1);
dim = size(qxy, 2);

dist = nan(n, 1);
vn = nan(n, dim);
idx_t = nan(n, 2);
for i = 1:n
    [tmp_dist, tmp_vn, tmp_t] = p2line_dist(qxy(i, :), xy(1:end-1, :), xy(2:end, :));
    [dist(i), tmp_idx] = min(tmp_dist);
    vn(i, :) = tmp_vn(tmp_idx, :);
    idx_t(i, :) = [tmp_idx, tmp_t(tmp_idx)];
end
end


function [dist, vp_n, t] = p2line_dist(p, xy1, xy2)
% INPUT
%   p:          1*2
%   xy1, xy2:   n*2

vp = bsxfun(@minus, p, xy1);
v2 = xy2 - xy1;

t = min(max(sum(vp .* v2, 2) ./ sum(v2.^2, 2), 0), 1);
vp_proj = bsxfun(@times, t, v2);
vp_n = bsxfun(@minus, vp, vp_proj);
dist = sqrt(sum(vp_n.^2, 2));
end


function [input_dist, output_dist] = input_output_distance(qx, curr_x, curr_j)
lon_offset = [-360, 1; 360, 1; 0, 1];
lat_offset = [-180, -1; 180, -1; 0, 1];
roll_offset = [-360, 1; 360, 1; 0, 1];
comb_offset = zeros(27, 6);  % offset, scale
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

q_num = size(qx, 1);
input_dist = inf(q_num, 1);
output_dist = inf(q_num, 1);

for comb_i = 1:27
    tmp_qx = bsxfun(@plus, bsxfun(@times, qx, comb_offset(comb_i, 4:6)), comb_offset(comb_i, 1:3));
    [dist, vn, it] = distance_to_poly_line(tmp_qx, curr_x);
    input_dist = min(input_dist, dist);
    for i = 1:size(vn, 1)
        if it(i, 2) > 1e-4 && it(i, 2) < 1 - 1e-4
            tmp_dy = vn(i, :) * curr_j(:, :, it(i, 1))';
            output_dist(i) = min(norm(tmp_dy), output_dist(i));
        end
    end
end
end