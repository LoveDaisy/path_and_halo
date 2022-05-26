function test_polygon2d_intersection()
test_cases = {@suite1};
num = length(test_cases);

fprintf('Start testing for polygon2d_intersection...\n');
for i = 1:num
    fprintf('Testing suite %d/%d...\n', i, num);
    test_cases{i}();
    fprintf('suite %d/%d passed!\n', i, num);
end
end

% ================================================================================
function suite1()
vtx2 = [0, 0; 1, 0; 1, 1; 0, 1];

fprintf(' case 1 ... ');
vtx1 = [-1, .5; -1, - .5; .5, - .5; .5, .5];
vtx0 = [0, 0; .5, 0; .5, .5; 0, .5];

vtx = geo.polygon2d_intersection(vtx1, vtx2);
assert(all(abs(vtx(:) - vtx0(:)) < 1e-8));
fprintf('passed!\n');

fprintf(' case 2 ... ');
vtx1 = [0, -.7; 1.7, 1; 2, 0];
vtx0 = [.7, 0; 1, .3; 1, 0];

vtx = geo.polygon2d_intersection(vtx1, vtx2);
assert(all(abs(vtx(:) - vtx0(:)) < 1e-8));
fprintf('passed!\n');
end
