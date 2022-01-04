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
crystal = opt.make_prism_crystal(1);
trace.fid = [3; 5];

quat = [1, 0, 0, 0];
ray_in_ll = [-40 + 180, 0];
ray_out_ll0 = [161.846476470779, 0];

[ray_out_ll, g_rot] = opt.crystal_system(quat, ray_in_ll, crystal, trace);
assert(all(abs(ray_out_ll(:) - ray_out_ll0(:)) < 1e-8));

dq = 1e-8;
ray_out_ll1 = opt.crystal_system(quat + [dq, 0, 0, 0], ray_in_ll, crystal, trace);
ray_out_ll2 = opt.crystal_system(quat + [0, dq, 0, 0], ray_in_ll, crystal, trace);
ray_out_ll3 = opt.crystal_system(quat + [0, 0, dq, 0], ray_in_ll, crystal, trace);
ray_out_ll4 = opt.crystal_system(quat + [0, 0, 0, dq], ray_in_ll, crystal, trace);
g_rot0 = [ray_out_ll1 - ray_out_ll; ray_out_ll2 - ray_out_ll; ray_out_ll3 - ray_out_ll;
    ray_out_ll4 - ray_out_ll]' / dq;

g_rot = g_rot(1:end-1, :);
assert(all(abs(g_rot(:) - g_rot0(:)) < 1e-5));
end
