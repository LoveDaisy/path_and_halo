clear; close all; clc;

n_side = 2^5;
n_pix = nSide2nPix(n_side);

input_phi = nan(n_pix, 1);

for i = 1:n_pix
    curr_ang = pix2ang(n_side, i);
    input_phi(i) = curr_ang{1}(2);
end

figure(1); clf;
imagesc(reshape(input_phi, [4*n_side, 3*n_side])');
axis equal; axis tight;