function y = halo_vis_fun(x, lim)
% Convert linear intensity x into visual intensity y

a = min(lim);
b = max(lim);
y = log10(x * b ./ (x + b) + a);
end