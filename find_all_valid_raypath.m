clear; close all; clc;

crystal = opt.make_prism_crystal(.2);
face_num = length(crystal.face);

max_hits = 5;
possible_raypaths = -ones(face_num^max_hits, max_hits);

tic;
% Find all possible raypaths
for i = 1:face_num^max_hits
    s = dec2base(i-1, face_num, max_hits);
    for j = 1:length(s)
        fn = s(j) - '0' + 1;
        if j > 1 && fn == possible_raypaths(i, j-1)
            break;
        else
            possible_raypaths(i, j) = fn;
        end
    end
end
fprintf('All possible raypaths: %d\n', size(possible_raypaths, 1));

% Reduce to normal form
for i = 1:size(possible_raypaths, 1)
    curr_raypath = possible_raypaths(i, :);
    idx = curr_raypath > 0;
    possible_raypaths(i, idx) = opt.reduce_raypath(curr_raypath(idx), 3);
end
possible_raypaths = unique(possible_raypaths, 'rows');
fprintf('After reduced: %d\n', size(possible_raypaths, 1));

% Check if is valid
raypath_valid_idx = false(size(possible_raypaths, 1), 1);
parfor i = 1:size(possible_raypaths, 1)
    fprintf('  checking raypath %d/%d...\n', i, size(possible_raypaths, 1));
    curr_raypath = possible_raypaths(i, :);
    idx = curr_raypath > 0;
    raypath_valid_idx(i) = opt.check_raypath(crystal, curr_raypath(idx));
end
possible_raypaths = possible_raypaths(raypath_valid_idx, :);
fprintf('Total valid: %d\n', size(possible_raypaths, 1));
toc;