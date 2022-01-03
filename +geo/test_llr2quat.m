function test_llr2quat()
test_cases = {@suite1};
num = length(test_cases);

fprintf('Start testing for llr2quat...\n');
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
mat = rotz(90 + lon) * rotx(90 - lat) * rotz(roll);

v = [3, 4, 1.2];
v_new0 = v * mat';

[q, jac] = geo.llr2quat([lon, lat, roll]);
v_new = quatrotate(q, v);

assert(all(abs(v_new - v_new0) < 1e-6));

d = 1e-4;
[q1, ~] = geo.llr2quat([lon, lat, roll] + [d, 0, 0]);
[q2, ~] = geo.llr2quat([lon, lat, roll] + [0, d, 0]);
[q3, ~] = geo.llr2quat([lon, lat, roll] + [0, 0, d]);

jac0 = [q1 - q; q2 - q; q3 - q]' / d;
assert(all(abs(jac0(:) - jac(:)) < 1e-8));
end
