function test_make_prism_crystal()
test_cases = {@suite1};
num = length(test_cases);

fprintf('Start testing for make_prism_crystal...\n');
for i = 1:num
    fprintf('Testing suite %d/%d...\n', i, num);
    test_cases{i}();
    fprintf('suite %d/%d passed!\n', i, num);
end
end

% ================================================================================
function suite1()
h = 1;
vtx0 = zeros(12, 3);
for i = 1:6
    vtx0(i, :) = [cosd(i * 60 - 90), sind(i * 60 - 90), h] / 2;
end
for i = 7:12
    vtx0(i, :) = [cosd(i * 60 - 90), sind(i * 60 - 90), -h] / 2;
end

face_norm0 = [0, 0, 1;
        0, 0, -1;
        1, 0, 0;
        1/2, sqrt(3) / 2, 0;
        - 1/2, sqrt(3) / 2, 0;
        - 1, 0, 0;
        - 1/2, -sqrt(3) / 2, 0;
        1/2, -sqrt(3) / 2, 0];

face_area0 = [3 * sqrt(3) / 8;
        3 * sqrt(3) / 8;
        h / 2;
        h / 2;
        h / 2;
        h / 2;
        h / 2;
        h / 2];

crystal = opt.make_prism_crystal(h);
assert(abs(h - crystal.h) < 1e-8);
assert(all(abs(vtx0(:) - crystal.vtx(:)) < 1e-8));
assert(all(abs(face_norm0(:) - crystal.face_norm(:)) < 1e-8));
assert(all(abs(face_area0(:) - crystal.face_area(:)) < 1e-8));
end
