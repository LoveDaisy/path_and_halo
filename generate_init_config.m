function config = generate_init_config(face_norm, n, level)
face_num = size(face_norm, 1);
n = [1; n(:)];

% initial grid
[r0, r0_ll, dr] = generate_healpix_grids(level);
valid_idx = r0 * face_norm(1, :)' < 0;

r1 = r0;
for i = 1:face_num
    r1(valid_idx, :) = refract_with_gradient(r1(valid_idx, :), face_norm(i, :), n(i), n(i+1));
    if i > 1
        valid_idx = valid_idx & (r1 * face_norm(i, :)' > 0);
    end
end
valid_idx = valid_idx & sum(r1.^2, 2) > 1e-4;
r1(~valid_idx, :) = nan;

bending_angle = acosd(sum(r0 .* r1, 2));
bending_angle_max = max(bending_angle);
bending_angle_min = min(bending_angle);

config.bending_angle = bending_angle;
config.bending_angle_max = bending_angle_max;
config.bending_angle_min = bending_angle_min;
config.r0 = r0;
config.r0_ll = r0_ll;
config.dr = dr;
end