function halo_img = generate_halo_image(lon_range, lat_range, res)
% Generate necessary data for halo image

halo_img_x = lon_range(1):res:lon_range(2);
halo_img_y = lat_range(1):res:lat_range(2);

img = nan(length(halo_img_y), length(halo_img_x));

halo_img.img = img;
halo_img.res = res;
halo_img.img_x = halo_img_x;
halo_img.img_y = halo_img_y;
halo_img.x_length = length(halo_img_x);
halo_img.y_length = length(halo_img_y);
end
