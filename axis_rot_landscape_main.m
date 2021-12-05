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

crystal.face_norm = face_norm;
trace.fid = [3; 5];
trace.n = [n; 1];

% initial grid
[axis_u, axis_u_ll, dr] = generate_healpix_grids(3);
u_num = size(axis_u, 1);
roll = (0:(2*dr):360)';
roll_num = length(roll);
total_num = u_num * roll_num;

sun_altitude = 20;
sun_longitude = 0;
ray_in = ll2xyz([sun_longitude, sun_altitude]);
ray_out_ll = zeros(total_num, 2);
jacob_store = zeros(2, 3, total_num);
axis_rot = nan(total_num, 3);

progress_cnt = 0;
progress_bin = 0.005;
for i = 1:total_num
    progress_cnt = progress_cnt + 1 / total_num;
    if progress_cnt > progress_bin
        progress_cnt = progress_cnt - floor(progress_cnt / progress_bin) * progress_bin;
        fprintf('processing %05.2f%%...\n', i / total_num * 100);
    end
    ll_idx = floor((i - 1) / roll_num) + 1;
    roll_idx = i - (ll_idx - 1) * roll_num;
    curr_rot = [axis_u_ll(ll_idx, :), roll(roll_idx)];
    axis_rot(i, :) = curr_rot;
    
    [tmp_ll, tmp_jcb] = crystal_system(curr_rot, [sun_longitude, sun_altitude], ...
        crystal, trace);
    ray_out_ll(i, :) = tmp_ll;
    jacob_store(:, :, i) = tmp_jcb;
end
ray_out = ll2xyz(ray_out_ll);
bending_angle = acosd(ray_out * ray_in');

%%
xx = kron(axis_u_ll(:, 1), ones(size(roll)));
yy = kron(axis_u_ll(:, 2), ones(size(roll)));
zz = repmat(roll, [u_num, 1]);

target_ll = [0, -5];
dd_ll = bsxfun(@minus, ray_out_ll, target_ll);
valid_idx = bending_angle > 0;
valid_idx = valid_idx & max(abs(dd_ll), [], 2) < dr;

uu = zeros(total_num, 1);
vv = zeros(total_num, 1);
ww = zeros(total_num, 1);
for i = 1:total_num
    tmp_uvw = dd_ll(i, :) * jacob_store(:, :, i);
    uu(i) = tmp_uvw(1);
    vv(i) = tmp_uvw(2);
    ww(i) = tmp_uvw(3);
end

figure(1); clf;
hold on;
scatter3(zz(valid_idx), xx(valid_idx), yy(valid_idx), 20, ...
    bending_angle(valid_idx), 'fill');
quiver3(zz(valid_idx), xx(valid_idx), yy(valid_idx), ...
    ww(valid_idx), uu(valid_idx), vv(valid_idx));
axis equal;
set(gca, 'xlim', [0, 360], 'ylim', [-180, 180], 'zlim', [-90, 90]);
colorbar;


figure(2); clf;
for i = 1:roll_num
    valid_idx = i:roll_num:total_num;
    scatter3(zz(valid_idx), xx(valid_idx), yy(valid_idx), 20, ...
        log10(sum(bsxfun(@minus, ray_out_ll(valid_idx, :), [0, -5]).^2, 2)), 'fill');
    set(gca, 'xlim', [0, 360], 'ylim', [-180, 180], 'zlim', [-90, 90]);
    drawnow;
    pause(.1);
end
