function [xyz, lon_lat, dr] = generate_healpix_grids(n_level)
% INPUT
% 	n_level:		n level for side, n_side = 2^n_level
% OUTPUT
%		xyz:				m*3, [x, y, z] unit vector on sphere
%		lon_lat:		m*2, [longitude, latitude] in degree
%		dr:					approximate angular resolution between two grid points. in degree

n_side = 2^n_level;
n_pix = nSide2nPix(n_side);
dr = sqrt(4 * pi / n_pix) * 180 / pi * 0.5;

[x, y, z] = pix2vec(n_side, 1:n_pix);
xyz = [x', y', z'];
lon_lat = [atan2d(y(:), x(:)), asind(z(:) ./ sqrt(sum(xyz.^2, 2)))];
end
