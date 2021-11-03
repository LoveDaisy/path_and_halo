function [v, g_rot] = rotation_with_gradient(v, lon_lat_theta)
% INPUT
%   v:              n*3, input vector
%   lon_lat_theta:  3-vector, [longitude, latitude, theta]

[u, g_u] = ll2xyz_with_gradient([lon_lat_theta(1), lon_lat_theta(2)]);

ux = [0, -u(3), u(2);
    u(3), 0, -u(1);
    -u(2), u(1), 0];
g_ux1 = [0, -g_u(3, 1), g_u(2, 1);
    g_u(3, 1), 0, -g_u(1, 1);
    -g_u(2, 1), g_u(1, 1), 0];
g_ux2 = [0, -g_u(3, 2), g_u(2, 2);
    g_u(3, 2), 0, -g_u(1, 2);
    -g_u(2, 2), g_u(1, 2), 0];

utu = u' * u;
g_utu1 = g_u(:, 1) * u + u' * g_u(:, 1)';
g_utu2 = g_u(:, 2) * u + u' * g_u(:, 2)';

cos_q = cosd(lon_lat_theta(3));
sin_q = sind(lon_lat_theta(3));
matRt = cos_q * eye(3) + sin_q * ux' + (1 - cos_q) * utu;
g_matRt1 = sin_q * g_ux1' + (1 - cos_q) * g_utu1;
g_matRt2 = sin_q * g_ux2' + (1 - cos_q) * g_utu2;
g_matRt3 = (-sin_q * eye(3) + cos_q * ux' + sin_q * utu) * pi / 180;

g_rot = [v * g_matRt1; v * g_matRt2; v * g_matRt3]';
v = v * matRt;
end