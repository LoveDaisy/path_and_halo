function test_distance_to_polyline()
test_cases = {@suite1};
num = length(test_cases);

fprintf('Start testing for distance_to_polyline...\n');
for i = 1:num
    fprintf('Testing suite %d/%d...\n', i, num);
    test_cases{i}();
    fprintf('suite %d/%d passed!\n', i, num);
end
end

function suite1()
line_pts = [0, 0;
        1, 0;
        1, 1;
        0, 1];

% Case 1
fprintf(' case 1...');
query_p = [0.5, 0.1];
d = 0.1;
v =- [0, 0.1];
it = [1, 0.5];
[tmp_d, tmp_v, tmp_it] = geo.distance_to_polyline(line_pts, query_p);

assert(all(abs(tmp_d - d) < 1e-8));
assert(all(abs(tmp_v - v) < 1e-8));
assert(all(abs(tmp_it - it) < 1e-8));
fprintf(' passed!\n');

% Case 2
fprintf(' case 2...');
query_p = [-0.5, 0.1];
d = norm([-0.5, .1]);
v = [0.5, -0.1];
it = [1, 0];
[tmp_d, tmp_v, tmp_it] = geo.distance_to_polyline(line_pts, query_p);

assert(all(abs(tmp_d - d) < 1e-8));
assert(all(abs(tmp_v - v) < 1e-8));
assert(all(abs(tmp_it - it) < 1e-8));
fprintf(' passed!\n');
end
