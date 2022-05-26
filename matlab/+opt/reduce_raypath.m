function raypath = reduce_raypath(raypath, sym_level)
% Reduce a raypath to a normal form. A normal form of a raypath is the one which has the smallest
% lexicographical order, up to certain level of symmetry.
% Applicable symmetries are (from lower to higher), none, C6, D6, D6h
%
% INPUT
%   raypath:        n*1
%   sym_level:      integer, must be one of 0, 1, 2, 3

p_idx = find(raypath > 2);
b_idx = find(raypath <= 2);

% None symmetry
if sym_level == 0
    return;
end

% Apply C6 symmetry
if ~isempty(p_idx)
    raypath(p_idx) = raypath(p_idx) - 3;
    raypath(p_idx) = mod(raypath(p_idx) - raypath(p_idx(1)), 6);
    raypath(p_idx) = raypath(p_idx) + 3;
end
if sym_level == 1
    return;
end

% Apply D6 symmetry
if ~isempty(p_idx)
    raypath(p_idx) = raypath(p_idx) - 3;
    tmp_path_pos = mod(raypath(p_idx), 6);
    tmp_path_neg = mod(-raypath(p_idx), 6);
    for i = 1:length(tmp_path_pos)
        if tmp_path_pos(i) < tmp_path_neg(i)
            raypath(p_idx) = tmp_path_pos;
            break;
        elseif tmp_path_pos(i) > tmp_path_neg(i)
            raypath(p_idx) = tmp_path_neg;
            break;
        end
    end
    raypath(p_idx) = raypath(p_idx) + 3;
end
if sym_level == 2
    return;
end

% Apply D6h symmetry
if ~isempty(b_idx)
    raypath(b_idx) = raypath(b_idx) - 1;
    raypath(b_idx) = mod(raypath(b_idx) - raypath(b_idx(1)), 2);
    raypath(b_idx) = raypath(b_idx) + 1;
end
end