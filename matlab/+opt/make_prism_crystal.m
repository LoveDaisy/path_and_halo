function crystal = make_prism_crystal(h)
% Make a prism crystal.
% The prism face 3 is towards x-axis, and the first vertex locates at -30 degree angle.
%
% INPUT
%   h:          scalar, a relative height for crystal. If the diameter is 1, the height is h.
%
% OUTPUT
%   crystal:    struct

crystal.n = 1.31;
crystal.h = h;

crystal.vtx = zeros(12, 3);
for i = 1:6
    crystal.vtx(i, :) = [cosd(i * 60 - 90), sind(i * 60 - 90), crystal.h] / 2;
end
for i = 7:12
    crystal.vtx(i, :) = [cosd(i * 60 - 90), sind(i * 60 - 90), -crystal.h] / 2;
end

crystal.face = {[1; 2; 3; 4; 5; 6];
            [12; 11; 10; 9; 8; 7];
            [1; 7; 8; 2];
            [2; 8; 9; 3];
            [3; 9; 10; 4];
            [4; 10; 11; 5];
            [5; 11; 12; 6];
            [6; 12; 7; 1]};
face_cnt = length(crystal.face);

crystal.face_norm = zeros(face_cnt, 3);
for i = 1:face_cnt
    tmp_v1 = diff(crystal.vtx(crystal.face{i}(1:2), :));
    tmp_v2 = diff(crystal.vtx(crystal.face{i}(2:3), :));
    tmp_norm = cross(tmp_v1, tmp_v2);
    tmp_norm = tmp_norm / norm(tmp_norm);
    crystal.face_norm(i, :) = tmp_norm;
end

crystal.face_area = zeros(face_cnt, 1);
for i = 1:face_cnt
    tmp_vtx = crystal.vtx(crystal.face{i}, :);
    tmp_vtx_cnt = size(tmp_vtx, 1);
    tmp_area = 0;
    for j = 2:tmp_vtx_cnt - 1
        tmp_v1 = tmp_vtx(j, :) - tmp_vtx(1, :);
        tmp_v2 = tmp_vtx(j + 1, :) - tmp_vtx(1, :);
        tmp_area = tmp_area + norm(cross(tmp_v1, tmp_v2)) / 2;
    end
    crystal.face_area(i) = tmp_area;
end
end
