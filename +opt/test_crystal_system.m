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
g_rot0 = [-0.0156722514735974, 0, -0.0156722514735974;
    0, 0.331223380243738, 0];

[ray_out_ll, g_rot] = opt.crystal_system(rot_llr, ray_in_ll, crystal, trace);
assert(all(abs(ray_out_ll(:) - ray_out_ll0(:)) < 1e-8));

assert(all(abs(g_rot(:) - g_rot0(:)) < 1e-8));
fprintf('passed!\n');

% ----------
fprintf(' case 2 ... ');
quat = [1, 0, 0, 0];
ray_in_ll = [-40 + 180, 0];
ray_out_ll0 = [161.846476470779, 0, 1];
g_rot0 = [0, 0, 0, 1.79590772980963;
    0, 37.9554035280461, -21.1055480272688, 0;
    1, 0, 0, 0];

[ray_out_ll, g_rot] = opt.crystal_system(quat, ray_in_ll, crystal, trace);
assert(all(abs(ray_out_ll(:) - ray_out_ll0(:)) < 1e-8));

assert(all(abs(g_rot(:) - g_rot0(:)) < 1e-8));
fprintf('passed!\n');
end

% ================================================================================
function suite2()
crystal = opt.make_prism_crystal(1);
trace.fid = [1; 3; 2; 4; 5; 1];

quat = [0.424935669519169, 0.480586073753202, -0.669476669301674, -0.374523285804727];
ray_in_ll = [180, -15];
ray_out_ll0 = [5, 25, 1];
g_rot0 = [39.0226918875199, 76.6161096090852, 48.9316069885360, 55.1212715447521;
    87.1735640313968, 132.680402326839, 115.713310660609, 62.3196193241224;
    0.424935669547968, 0.480586073785773, -0.669476669347046, -0.374523285830109];

[ray_out_ll, g_rot] = opt.crystal_system(quat, ray_in_ll, crystal, trace);
assert(all(abs(ray_out_ll(:) - ray_out_ll0(:)) < 1e-8));

assert(all(abs(g_rot(:) - g_rot0(:)) < 1e-8));
end
