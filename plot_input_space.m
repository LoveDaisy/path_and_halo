function plot_input_space(config, target_diff, seeds_idx, checked_idx, curr_x, curr_j)
hold on;
axis_rot = config.axis_rot_store;
curr_j = curr_j * 5;
scatter3(axis_rot(seeds_idx, 1), axis_rot(seeds_idx, 2), axis_rot(seeds_idx,3), 50, ...
    target_diff(seeds_idx), 'fill');
if ~isempty(curr_x)
    plot3(curr_x(:, 1), curr_x(:, 2), curr_x(:,3), 'r-s');
    for i = 1:size(curr_j, 3)
        plot3([0, curr_j(1, 1, i)] + curr_x(i, 1), [0, curr_j(1, 2, i)] + curr_x(i, 2), ...
            [0, curr_j(1, 3, i)] + curr_x(i, 3), 'm');
        plot3([0, curr_j(2, 1, i)] + curr_x(i, 1), [0, curr_j(2, 2, i)] + curr_x(i, 2), ...
            [0, curr_j(2, 3, i)] + curr_x(i, 3), 'm');
        
        m = curr_j(:, :, i) * curr_j(:, :, i)';
        m2 = m * m;
        cos_q = cosd(0:10:360)';
        sin_q = sind(0:10:360)';
        v2 = sqrt(m2(1,1)/(m2(1,1)*m2(2,2) - m2(1,2)^2)) * sin_q;
        v1 = cos_q / sqrt(m2(1,1)) - m2(1,2)/m2(1,1) * v2;
        u = [v1, v2] * curr_j(:, :, i);
        u = bsxfun(@plus, u, curr_x(i, :));
        plot3(u(:, 1), u(:, 2), u(:, 3), 'k');
    end
end
plot3(axis_rot(checked_idx, 1), axis_rot(checked_idx, 2), axis_rot(checked_idx,3), 'cx');
axis equal;
colorbar;
end