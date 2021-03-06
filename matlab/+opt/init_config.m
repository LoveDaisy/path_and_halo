function config = init_config(crystal, trace, sun_ll, varargin)
% Initialize a configuration for calculation
%
% INPUT
%   crystal:        struct
%   trace:          struct
%   sun_ll:         [longitude, latitude]
%
% OPTIONS
%   'GridLevel':    integer, default is 3, the level for HEALPix generation

p = inputParser;
p.addParameter('GridLevel', 3, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
p.parse(varargin{:});

% initial grid
[~, rot0_ll, dr] = geo.generate_healpix_grids(p.Results.GridLevel);
roll0 = (0:dr * 2:360)';

a_num = size(rot0_ll, 1);
r_num = length(roll0);
total_num = a_num * r_num;

axis_llr_store = [repmat(rot0_ll, [r_num, 1]), kron(roll0, ones(a_num, 1))];

out_xyz = zeros(total_num, 3);
out_ll = zeros(total_num, 2);

progress_cnt = 0;
progress_bin = 0.1;
ray_in_ll = [sun_ll(1) + 180, -sun_ll(2)];
for i = 1:total_num
    if progress_cnt > progress_bin
        fprintf('init_config (%04.1f%%)...\n', i / total_num * 100);
        progress_cnt = progress_cnt - floor(progress_cnt / progress_bin) * progress_bin;
    end
    progress_cnt = progress_cnt + 1 / total_num;

    tmp_out_xyz = opt.crystal_system(axis_llr_store(i, :), ray_in_ll, crystal, trace);
    out_xyz(i, :) = tmp_out_xyz;
    out_ll(i, :) = geo.xyz2ll(tmp_out_xyz);
end

in_xyz = geo.ll2xyz(ray_in_ll);
bending_angle = acosd(out_xyz * in_xyz');

config.bending_angle = bending_angle;
config.axis_llr_store = axis_llr_store;
config.axis_quat_store = geo.llr2quat(axis_llr_store);
config.dr = dr;
config.out_ll = out_ll;
config.out_xyz = out_xyz;
config.sun_ll = sun_ll;
config.ray_in_ll = ray_in_ll;
config.crystal = crystal;
config.trace = trace;
end
