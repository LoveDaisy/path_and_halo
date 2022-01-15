function show_halo_img(halo_img)
% Display halo image.
%
% INPUT
%   halo_img:       struct, generated from generate_halo_img()

vis_range = prctile(halo_img.img(:), [10, 90]);
vis_range = max(vis_range, 1e-3 * vis_range(2));
imagesc(halo_img.img_x, halo_img.img_y, halo_vis_fun(halo_img.img, vis_range));
axis equal; axis tight; axis xy;
end

function y = halo_vis_fun(x, lim)
% Convert linear intensity x into visual intensity y

a = min(lim);
b = max(lim);
y = log10(x * b ./ (x + b) + a);
end
