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
crystal = opt.make_prism_crystal(1);
trace.fid = [3; 5];

% ----------
fprintf(' case 1 ... ');
rot_llr = [-90, 90, 0];
ray_in_ll = [-40 + 180, 0];
ray_out_ll0 = [161.846476470779, 0];
ray_out_xyz0 = geo.ll2xyz(ray_out_ll0);
g_rot0 = [8.52182180377525e-05,2.17603712826531e-09,8.52182191479756e-05;
    0.000259903028543285,1.27453603226968e-09,0.000259903027433062;
    0,0.00578093854503393,0];

[ray_out_xyz, g_rot] = opt.crystal_system(rot_llr, ray_in_ll, crystal, trace);
assert(all(abs(ray_out_xyz(:) - ray_out_xyz0(:)) < 1e-8));

assert(all(abs(g_rot(:) - g_rot0(:)) < 3e-7));
fprintf('passed!\n');

% ----------
fprintf(' case 2 ... ');
quat = [1, 0, 0, 0];
ray_in_ll = [-40 + 180, 0];
ray_out_ll0 = [161.846476470779, 0];
ray_out_xyz0 = [geo.ll2xyz(ray_out_ll0), 1];
g_rot0 = [0,0,0,-0.00976582590801034;
    0,0,0,-0.0297843332868331;
    0,0.662446760487477,-0.368361303513634,0;
    1,0,0,0];

[ray_out_xyz, g_rot] = opt.crystal_system(quat, ray_in_ll, crystal, trace);
assert(all(abs(ray_out_xyz(:) - ray_out_xyz0(:)) < 1e-8));

assert(all(abs(g_rot(:) - g_rot0(:)) < 1e-8));
fprintf('passed!\n');
end

% ================================================================================
function suite2()
crystal = opt.make_prism_crystal(1);
trace.fid = [1; 3; 2; 4; 5; 1];

quat = [0.424935669519169, 0.480586073753202, -0.669476669301674, -0.374523285804727];
ray_in_ll = [180, -15];
ray_out_ll0 = [5, 25];
ray_out_xyz0 = [geo.ll2xyz(ray_out_ll0), 1];
g_rot0 = [-0.694350411387566,-1.08056277836721,-0.917721580714424,-0.533917383562988;
    0.558873138601976,1.12201016546215,0.696669177227574,0.828530141768223;
    1.37891622344290,2.09874589084732,1.83035950303392,0.985775161078285;
    0.424935669547968,0.480586073785773,-0.669476669347046,-0.374523285830109];

[ray_out_xyz, g_rot] = opt.crystal_system(quat, ray_in_ll, crystal, trace);
assert(all(abs(ray_out_xyz(:) - ray_out_xyz0(:)) < 1e-8));

assert(all(abs(g_rot(:) - g_rot0(:)) < 1e-8));
end
