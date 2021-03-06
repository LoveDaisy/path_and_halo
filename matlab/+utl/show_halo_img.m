function show_halo_img(halo_img, dr, w)
% Display halo image.
%
% INPUT
%   halo_img:       struct, generated from generate_halo_image()

if nargin == 1
    dr = 4;
    w = 95;
elseif nargin == 2
    w = 95;
end
k = 2^(-dr);

white_lim = prctile(halo_img.img(:), w);
black_lim = k * white_lim;
if dr < 0
    imagesc(halo_img.img_x, halo_img.img_y, halo_img.img);
else
    imagesc(halo_img.img_x, halo_img.img_y, vis_fun(halo_img.img, [black_lim, white_lim]));
end
axis equal; axis tight; axis xy;
end

function y = vis_fun(x, lim)
y = (x + lim(1)) * lim(2) ./ (x + lim(2));
y = log10(y);
end
