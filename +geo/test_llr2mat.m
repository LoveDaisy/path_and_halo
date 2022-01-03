function test_llr2mat()
test_cases = {@suite1};
num = length(test_cases);

fprintf('Start testing for llr2mat...\n');
for i = 1:num
    fprintf('Testing suite %d/%d...\n', i, num);
    test_cases{i}();
    fprintf('suite %d/%d passed!\n', i, num);
end
end

% ================================================================================
function suite1()
lon = 10;
lat = 35;
roll = 60;
mat0 = rotz(90 + lon) * rotx(90 - lat) * rotz(roll);

[mat, jac] = geo.llr2mat([lon, lat, roll]);
assert(all(abs(mat(:) - mat0(:)) < 1e-10));

d = 1e-3;
mat1 = geo.llr2mat([lon, lat, roll] + [d, 0, 0]);
mat2 = geo.llr2mat([lon, lat, roll] + [0, d, 0]);
mat3 = geo.llr2mat([lon, lat, roll] + [0, 0, d]);
jac0 = cat(3, (mat1 - mat0) / d, (mat2 - mat0) / d, (mat3 - mat0) / d);
assert(all(abs(jac(:) - jac0(:)) < 5e-7));
end
