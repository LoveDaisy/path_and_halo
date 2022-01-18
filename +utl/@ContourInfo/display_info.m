function display_info(obj)
line_color = [0, 0.447, 0.741;
        0.85, 0.325, 0.098;
        0.929, 0.694, 0.125;
        0.494, 0.184, 0.556;
        0.466, 0.674, 0.188;
        0.301, 0.745, 0.933;
        0.635, 0.078, 0.184;
        0, 0.447, 0.741;
        0.85, 0.325, 0.098;
        0.929, 0.694, 0.125];

subplot(2, 2, 1);
hold on;
for k = 1:obj.total_cnt
    plot(obj.llr_contour_store{k}(:, 2), obj.llr_contour_store{k}(:, 3), '-o');
    plot(obj.llr_interp_store{k}(:, 2), obj.llr_interp_store{k}(:, 3), '.-');
end
box on;
roll_lim1 = get(gca, 'ylim');
subplot(2, 2, 2);
hold on;
for k = 1:obj.total_cnt
    plot(obj.llr_contour_store{k}(:, 1), obj.llr_contour_store{k}(:, 3), '-o');
    plot(obj.llr_interp_store{k}(:, 1), obj.llr_interp_store{k}(:, 3), '.-');
end
box on;
roll_lim2 = get(gca, 'ylim');
subplot(2, 2, 1);
set(gca, 'ylim', [min(roll_lim1(1), roll_lim2(1)), max(roll_lim1(2), roll_lim2(2))]);
subplot(2, 2, 2);
set(gca, 'ylim', [min(roll_lim1(1), roll_lim2(1)), max(roll_lim1(2), roll_lim2(2))]);

subplot(2, 2, [3, 4]);
hold on;
for k = 1:obj.total_cnt
    plot(obj.weight_component_store{k}(:, 1), obj.weight_component_store{k}(:, 2), '.-k', 'linewidth', 2);
    plot(obj.weight_component_store{k}(:, 1), obj.weight_component_store{k}(:, 3), '.-', ...
        'color', line_color(1, :));
    plot(obj.weight_component_store{k}(:, 1), obj.weight_component_store{k}(:, 4), '--', 'linewidth', 2, ...
        'color', line_color(2, :));
    plot(obj.weight_component_store{k}(:, 1), obj.weight_component_store{k}(:, 5), ':', 'linewidth', 3, ...
        'color', line_color(3, :));
    plot(obj.weight_component_store{k}(:, 1), obj.weight_component_store{k}(:, 6), '-', 'linewidth', 1.2, ...
        'color', line_color(4, :));
end
plot([obj.weight_component_store{1}(1, 1), obj.weight_component_store{1}(end, 1)], [1, 1], ':k');
set(gca, 'yscale', 'log', 'ylim', [1e-10, 1e4]);
box on;
end
