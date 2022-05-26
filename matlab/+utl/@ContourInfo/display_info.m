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

clf;
fig_pos = get(gcf, 'position');

outer_spacing = 0.05;
inner_spacing = 0.02;
x0 = outer_spacing;
y0 = outer_spacing;
h1 = 0.4; w1 = h1 * fig_pos(4) / fig_pos(3);
h2 = 0.18;
h3 = 1 - 3 * outer_spacing - h1 - h2;

if obj.total_cnt < 1
    return;
end
data_dim = size(obj.contour_store{1}, 2);
if data_dim > 3
    subspace_idx = [2, 3, 4];
end

% ----------------------------------
% Plot candidate seeds and rot contours
axes('position', [x0, 1 - y0 - h1 + inner_spacing, w1 - inner_spacing, h1 - inner_spacing]);
hold on;
if ~isempty(obj.candidate_seeds)
    utl.plot_data_3d(obj.candidate_seeds, subspace_idx, 'o', 'color', line_color(1, :));
end
for i = 1:obj.total_cnt
    utl.plot_data_3d(obj.contour_store{i}, subspace_idx, 'x-', 'color', line_color(2, :));
end
box on;
axis equal;

% ----------------------------------
% Plot longitude-roll plane
a12 = axes('position', [x0 + w1, 1 - y0 - h1, 1 - 2 * x0 - w1, h1]);
hold on;
for i = 1:obj.total_cnt
    plot(obj.llr_contour_store{i}(:, 1), obj.llr_contour_store{i}(:, 3), 'o');
    plot(obj.llr_interp_store{i}(:, 1), obj.llr_interp_store{i}(:, 3), '.-');
end
box on;
set(a12, 'YAxisLocation', 'right', 'XAxisLocation', 'top');

% ----------------------------------
% Plot roll-latitude plane
a21 = axes('position', [x0, y0 + h3 + outer_spacing, w1, h2]);
hold on;
for i = 1:obj.total_cnt
    plot(obj.llr_contour_store{i}(:, 3), obj.llr_contour_store{i}(:, 2), 'o');
    plot(obj.llr_interp_store{i}(:, 3), obj.llr_interp_store{i}(:, 2), '.-');
end
set(a21, 'xdir', 'reverse');
box on;

% ----------------------------------
% Plot longitude-latitude plane
a22 = axes('position', [x0 + w1, y0 + h3 + outer_spacing, 1 - 2 * x0 - w1, h2]);
hold on;
for i = 1:obj.total_cnt
    plot(obj.llr_contour_store{i}(:, 1), obj.llr_contour_store{i}(:, 2), 'o');
    plot(obj.llr_interp_store{i}(:, 1), obj.llr_interp_store{i}(:, 2), '.-');
end
box on;
set(a22, 'YAxisLocation', 'right');

lon_lim = get(a22, 'xlim');
if abs(diff(lon_lim)) < 5
    lon_lim = [-2.5, 2.5] + mean(lon_lim);
    set(a22, 'xlim', lon_lim);
    set(a12, 'xlim', lon_lim);
end

lat_lim = get(a22, 'ylim');
if abs(diff(lat_lim)) < 5
    lat_lim = [-2.5, 2.5] + mean(lat_lim);
    set(a22, 'ylim', lat_lim);
    set(a21, 'ylim', lat_lim);
end

roll_lim2 = get(a12, 'ylim');
roll_lim1 = get(a21, 'xlim');
roll_lim = [min(roll_lim1(1), roll_lim2(2)), max(roll_lim1(2), roll_lim2(2))];
set(a12, 'ylim', roll_lim);
set(a21, 'xlim', roll_lim);

% ----------------------------------
% Plot weight components
axes('position', [x0, y0, 1 - 2 * x0, h3]);
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
xlim = get(gca, 'xlim');
plot(xlim, [1, 1], ':k');
set(gca, 'yscale', 'log', 'ylim', [1e-8, 1e3]);
box on;
end
