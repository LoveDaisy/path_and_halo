function test_crystal_system()
test_cases = {@suite1, @suite2};
num = length(test_cases);

fprintf('Start testing for crystal_system...\n');
for i = 1:num
    fprintf('Testing suite %d/%d...\n', i, num);
    test_cases{i}();
    fprintf('suite %d/%d passed!\n', i, num);
end
end

% ================================================================================
function suite1()
crystal.face_norm = [0, 0, 1; % face 1
                0, 0, -1; % face 2
                1, 0, 0; % face 3
                1/2, sqrt(3) / 2, 0; % face 4
                -1/2, sqrt(3) / 2, 0; % face 5
                -1, 0, 0; % face 6
                -1/2, -sqrt(3) / 2, 0; % face 7
                1/2, -sqrt(3) / 2, 0]; % face 8
crystal.h = 1;
crystal.vtx = zeros(12, 3);
for i = 1:6
    crystal.vtx(i, :) = [cosd(i * 60 - 30), sind(i * 60 - 30), crystal.h];
end
for i = 7:12
    crystal.vtx(i, :) = [cosd(i * 60 - 30), sind(i * 60 - 30), -crystal.h];
end
crystal.face = {[6; 5; 4; 3; 2; 1];
            [7; 8; 9; 10; 11; 12];
            [1; 2; 8; 7];
            [2; 3; 9; 8];
            [3; 4; 10; 9];
            [4; 5; 11; 10];
            [5; 6; 12; 11];
            [6; 1; 7; 12]};
crystal.face_area = [3 * sqrt(3) / 2; 3 * sqrt(3) / 2; crystal.h * ones(6, 1) * 2];
crystal.n = 1.31;

trace.fid = [3; 5];

rot_llr = [-90, 90, 0];
ray_in_ll = [-40 + 180, 0];
ray_out_ll0 = [161.846476470779, 0];

[ray_out_ll, g_rot] = opt.crystal_system(rot_llr, ray_in_ll, crystal, trace);
assert(all(abs(ray_out_ll(:) - ray_out_ll0(:)) < 1e-8));

dq = 1e-3;
ray_out_ll1 = opt.crystal_system(rot_llr + [dq, 0, 0], ray_in_ll, crystal, trace);
ray_out_ll2 = opt.crystal_system(rot_llr + [0, dq, 0], ray_in_ll, crystal, trace);
ray_out_ll3 = opt.crystal_system(rot_llr + [0, 0, dq], ray_in_ll, crystal, trace);
g_rot0 = [ray_out_ll1 - ray_out_ll; ray_out_ll2 - ray_out_ll; ray_out_ll3 - ray_out_ll]' / dq;
assert(all(abs(g_rot(:) - g_rot0(:)) < 1e-5));
end

% ================================================================================
function suite2()
crystal.face_norm = [0, 0, 1; % face 1
                0, 0, -1; % face 2
                1, 0, 0; % face 3
                1/2, sqrt(3) / 2, 0; % face 4
                -1/2, sqrt(3) / 2, 0; % face 5
                -1, 0, 0; % face 6
                -1/2, -sqrt(3) / 2, 0; % face 7
                1/2, -sqrt(3) / 2, 0]; % face 8
crystal.h = 1;
crystal.vtx = zeros(12, 3);
for i = 1:6
    crystal.vtx(i, :) = [cosd(i * 60 - 30), sind(i * 60 - 30), crystal.h];
end
for i = 7:12
    crystal.vtx(i, :) = [cosd(i * 60 - 30), sind(i * 60 - 30), -crystal.h];
end
crystal.face = {[6; 5; 4; 3; 2; 1];
            [7; 8; 9; 10; 11; 12];
            [1; 2; 8; 7];
            [2; 3; 9; 8];
            [3; 4; 10; 9];
            [4; 5; 11; 10];
            [5; 6; 12; 11];
            [6; 1; 7; 12]};
crystal.face_area = [3 * sqrt(3) / 2; 3 * sqrt(3) / 2; crystal.h * ones(6, 1) * 2];
crystal.n = 1.31;

trace.fid = [3; 5];

quat3 = [0, 0, 0];
ray_in_ll = [-40 + 180, 0];
ray_out_ll0 = [161.846476470779, 0];

[ray_out_ll, g_rot] = opt.crystal_system(quat3, ray_in_ll, crystal, trace, 'RotSpace', 'quat3');
assert(all(abs(ray_out_ll(:) - ray_out_ll0(:)) < 1e-8));

dq = 1e-8;
ray_out_ll1 = opt.crystal_system(quat3 + [dq, 0, 0], ray_in_ll, crystal, trace, 'RotSpace', 'quat3');
ray_out_ll2 = opt.crystal_system(quat3 + [0, dq, 0], ray_in_ll, crystal, trace, 'RotSpace', 'quat3');
ray_out_ll3 = opt.crystal_system(quat3 + [0, 0, dq], ray_in_ll, crystal, trace, 'RotSpace', 'quat3');
g_rot0 = [ray_out_ll1 - ray_out_ll; ray_out_ll2 - ray_out_ll; ray_out_ll3 - ray_out_ll]' / dq;
assert(all(abs(g_rot(:) - g_rot0(:)) < 1e-5));
end