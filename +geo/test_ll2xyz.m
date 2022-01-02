function test_ll2xyz()
test_cases = {@suite1, @suite2};
num = length(test_cases);

fprintf('Start testing for ll2xyz...\n');
for i = 1:num
    fprintf('Testing suite %d/%d...\n', i, num);
    test_cases{i}();
    fprintf('suite %d/%d passed!\n', i, num);
end
end

% ================================================================================
function suite1()
lon = 10;
lat = 20;
dq = 1e-4;
ll = [lon, lat];
xyz0 = [cosd(lat) * cosd(lon), cosd(lat) * sind(lon), sind(lat)];

tmp_xyz1 = [cosd(lat) * cosd(lon + dq), cosd(lat) * sind(lon + dq), sind(lat)];
tmp_xyz2 = [cosd(lat + dq) * cosd(lon), cosd(lat + dq) * sind(lon), sind(lat + dq)];
jac0 = [tmp_xyz1 - xyz0; tmp_xyz2 - xyz0]' / dq;

[xyz, jac] = geo.ll2xyz(ll);
assert(all(abs(xyz(:) - xyz0(:)) < 1e-6));
assert(all(abs(jac(:) - jac0(:)) < 1e-6));
end

% ================================================================================
function suite2()
rng(1234);
num = 100;
lon = rand(num, 1) * 360;
lat = rand(num, 1) * 180 - 90;
dq = 1e-4;
ll = [lon, lat];
xyz0 = [cosd(lat) .* cosd(lon), cosd(lat) .* sind(lon), sind(lat)];
jac0 = zeros(3, 2, num);

for i = 1:num
    tmp_xyz1 = [cosd(lat(i)) * cosd(lon(i) + dq), cosd(lat(i)) * sind(lon(i) + dq), sind(lat(i))];
    tmp_xyz2 = [cosd(lat(i) + dq) * cosd(lon(i)), cosd(lat(i) + dq) * sind(lon(i)), sind(lat(i) + dq)];
    jac0(:, :, i) = [tmp_xyz1 - xyz0(i, :); tmp_xyz2 - xyz0(i, :)]' / dq;
end

[xyz, jac] = geo.ll2xyz(ll);
assert(all(abs(xyz(:) - xyz0(:)) < 1e-6));
assert(all(abs(jac(:) - jac0(:)) < 1e-6));
end
