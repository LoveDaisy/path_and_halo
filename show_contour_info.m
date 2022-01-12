function show_contour_info(contour_store, rot_store, cmp_store)
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
for k = 1:length(contour_store)
    plot(contour_store{k}(:, 2), contour_store{k}(:, 3), '-o');
    plot(rot_store{k}(:, 2), rot_store{k}(:, 3), '.-');
end
box on;
roll_lim1 = get(gca, 'ylim');
subplot(2, 2, 2);
hold on;
for k = 1:length(contour_store)
    plot(contour_store{k}(:, 1), contour_store{k}(:, 3), '-o');
    plot(rot_store{k}(:, 1), rot_store{k}(:, 3), '.-');
end
box on;
roll_lim2 = get(gca, 'ylim');
subplot(2, 2, 1);
set(gca, 'ylim', [min(roll_lim1(1), roll_lim2(1)), max(roll_lim1(2), roll_lim2(2))]);
subplot(2, 2, 2);
set(gca, 'ylim', [min(roll_lim1(1), roll_lim2(1)), max(roll_lim1(2), roll_lim2(2))]);
subplot(2, 2, [3, 4]);
hold on;
for k = 1:length(cmp_store)
    plot(cmp_store{k}(:, 1), cmp_store{k}(:, 2), '.-k', 'linewidth', 2);
    plot(cmp_store{k}(:, 1), cmp_store{k}(:, 3), '.-', ...
        'color', line_color(1, :));
    plot(cmp_store{k}(:, 1), cmp_store{k}(:, 4), '--', 'linewidth', 2, ...
        'color', line_color(2, :));
    plot(cmp_store{k}(:, 1), cmp_store{k}(:, 5), ':', 'linewidth', 3, ...
        'color', line_color(3, :));
    plot(cmp_store{k}(:, 1), cmp_store{k}(:, 6), '-', 'linewidth', 1.2, ...
        'color', line_color(4, :));
end
plot([cmp_store{1}(1, 1), cmp_store{1}(end, 1)], [1, 1], ':k');
set(gca, 'yscale', 'log', 'ylim', [1e-8, 1e4]);
box on;
end
