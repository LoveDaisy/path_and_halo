clear; close all; clc;

face_norm = [0, 0, 1;     % face 1
    0, 0, -1;             % face 2
    -sqrt(3)/2, 1/2, 0;   % face 3
    0, 1, 0;              % face 4
    sqrt(3)/2, 1/2, 0;    % face 5
    sqrt(3)/2, -1/2, 0;   % face 6
    0, -1, 0;             % face 7
    -sqrt(3)/2, -1/2, 0]; % face 8
n = 1.31;

[r0, r0_ll, dr] = generate_healpix_grids(5);

r1 = nan(size(r0));
r2 = nan(size(r0));

entry_face_idx = 3;
exit_face_idx = 5;

crystal_zenith = [90, 0.2];  % mean, std
tmp_x = linspace(-90, 90, 5000);
tmp_pdf = exp(-(90 - tmp_x - crystal_zenith(1)).^2 / 2 / crystal_zenith(2)^2) / crystal_zenith(2);
axis_pdf = @(lon, lat) (exp(-(90 - lat - crystal_zenith(1)).^2 / 2 / crystal_zenith(2)^2) / ...
    crystal_zenith(2) / (sum(tmp_pdf) * (tmp_x(2) - tmp_x(1)))) * ...
    (1 / 360);

%%
sun_altitude = 20;
sun_longitude = 0;
ray_in_xyz = [cosd(sun_altitude) * cosd(sun_longitude), ...
    cosd(sun_altitude) * sind(sun_longitude), ...
    sind(sun_altitude)];

config = generate_init_config(face_norm([entry_face_idx, exit_face_idx], :), [1.31, 1], 6);

% halo_img_x = linspace(-180, 180, 512);
% halo_img_y = linspace(-90, 90, 256);
halo_img_x = -10:.5:10;
halo_img_y = -30:.5:0;
halo_img = zeros(length(halo_img_y), length(halo_img_x));
checked_pix = 0;
progress_cnt = 0;
progress_bin = 0.001;
for w = 1:length(halo_img_x)
    for h = 1:length(halo_img_y)
% for w = 240:256
%     for h = 90:150
% for w = 100
%     for h = 240
        lon = halo_img_x(w);
        lat = halo_img_y(h);
        ray_out_xyz = ll2xyz_with_gradient([lon, lat]);
        curr_bending_angle = acosd(dot(ray_out_xyz, ray_in_xyz));
        curr_target_input_output = [ray_in_xyz, ray_out_xyz];
        
        checked_pix = checked_pix + 1;
        progress_cnt = progress_cnt + 1 / numel(halo_img);
        
        [x_contour, g_angle, y_val, jacobian] = find_bending_angle_contour(curr_bending_angle, ...
            face_norm([entry_face_idx, exit_face_idx], :), [1.31, 1], 'config', config);
        if isempty(x_contour)
            continue;
        end
        
        weight = 0;
%         if abs(lon + 46.14) < 7e-3 && abs(lat - 33.53) < 7e-3
%             fprintf('!!!\n');
%         end
        p_store = {};
        axis_store = {};
        for k = 1:length(x_contour)
            curr_x_ll = x_contour{k};
            curr_y_ll = y_val{k};
            curr_g_angle = g_angle{k};
            curr_g_a_norm = sqrt(sum(curr_g_angle.^2, 2));
            
            curr_x_xyz = ll2xyz_with_gradient(curr_x_ll);
            curr_y_xyz = ll2xyz_with_gradient(curr_y_ll);
            curr_origin_input_output = [curr_x_xyz, curr_y_xyz];
            
            curr_q = find_rotation(curr_origin_input_output, ...
                repmat(curr_target_input_output, [size(curr_x_ll, 1), 1]));
            curr_axis = quatrotate(curr_q, [0, 0, 1]);
            curr_axis_ll = xyz2ll_with_gradient(curr_axis);
            
            p = axis_pdf(curr_axis_ll(:, 1), curr_axis_ll(:, 2));
            
            % then interpolation
            for kk = 2:size(curr_axis_ll, 1)
                if max(p(kk-1:kk)) < 1e-10
                    continue;
                end
                ds = acosd(dot(curr_axis(kk, :), curr_axis(kk-1, :)));
                interp_num = max(floor(ds / 0.05) + 1, 2);
                tmp_ds = linspace(0, ds, interp_num)';
                tmp_q = quatinterp(repmat(curr_q(kk-1, :), [interp_num, 1]), ...
                    repmat(curr_q(kk, :), [interp_num, 1]), tmp_ds / ds);
                tmp_axis = quatrotate(tmp_q, [0, 0, 1]);
                tmp_axis_ll = xyz2ll_with_gradient(tmp_axis);

                tmp_p = axis_pdf(tmp_axis_ll(:, 1), tmp_axis_ll(:, 2));
                tmp_norm = interp1([0; ds], curr_g_a_norm([kk-1; kk]), tmp_ds);
                tmp_p = tmp_p ./ tmp_norm;
                weight = weight + nansum((tmp_p(1:end-1) + tmp_p(2:end)) / 2 .* diff(tmp_ds));
            end
            
            % p = axis_pdf(curr_axis_ll(:, 1), curr_axis_ll(:, 2));
            
            % p = p ./ max(curr_g_a_norm, 1e-3);
            % weight = weight + nansum((p(1:end-1) + p(2:end)) / 2 .* diff(interp_s));
            % p_store{k} = [interp_s, p, [0; cumsum(sqrt(sum(diff(curr_axis_ll).^2, 2)))]];
            axis_store{k} = curr_axis_ll;
        end
        halo_img(h, w) = weight;
        
        if weight < 1e-8
            continue;
        end
        
%         if progress_cnt > progress_bin
            fprintf('progressing %04.1f%%...\n', checked_pix / numel(halo_img) * 100);
            progress_cnt = progress_cnt - floor(progress_cnt / progress_bin) * progress_bin;
            figure(4); clf;
            tmp_f4_pos = get(gcf, 'position');
            imagesc(halo_img_x, halo_img_y, halo_img);
            axis equal; axis tight; axis xy;
            drawnow;
            
            figure(3); clf;
            set(gcf, 'position', tmp_f4_pos + [tmp_f4_pos(3), 0, 0, 0]);
            hold on;
            for k = 1:length(x_contour)
                plot(x_contour{k}(:,1), x_contour{k}(:,2), '-o');
                plot(axis_store{k}(:, 1), axis_store{k}(:, 2), '-s');
            end
            title(sprintf('lon lat: (%d,%d) = (%05.2f,%05.2f)\nw: %05.2e', w, h, lon, lat, weight));
            axis equal; axis tight;
            box on;
            set(gca, 'xlim', [-180, 180], 'ylim', [-90, 90]);
            drawnow;
            
            % figure(6); clf;
            % set(gcf, 'position', tmp_f4_pos + [0, -tmp_f4_pos(4), 0, 0]);
            % hold on;
            % for k = 1:length(x_contour)
                % plot(p_store{k}(:, 1), p_store{k}(:, 2), '-o');
                % plot(p_store{k}(:, 3), p_store{k}(:, 2), '-s');
            % end
            % title(sprintf('w: %05.3e', weight));
            % box on;
            % set(gca, 'yscale', 'log', 'ylim', [1e-12, 1e-1]);
            % drawnow;
%         end
    end
end

%%
figure(4); clf;
imagesc(halo_img_x, halo_img_y, halo_img);
axis equal; axis tight; axis xy;
