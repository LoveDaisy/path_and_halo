function vtx = polygon2d_intersection(vtx1, vtx2)
% Calculate intersection of two *CONVEX* 2D polygon.
% If any of vtx1 or vtx2 is not convex, the result may not be right.
%
% INPUT
%   vtx1:       n*2, polygon vertexes
%   vtx2:       m*2, polygon vertexes
%
% OUTPUT
%   vtx:        k*2, intersection polygon vertexes

num1 = size(vtx1, 1);
num2 = size(vtx2, 1);
if num1 < 3 || num2 < 3
    error('polygon vertexes must be more than three!');
end

% Iterate every edge of polygon 2, and cut off outer part of polygon 1.
vtx = nan(max(num1, num2) * 2, 2);

% First find out which side is inner side. Those have same sign to z0 considered as inner.
z0 = det(diff(vtx2(1:3, :)));

eps_z = 1e-8;
% Then start cutting
for i = 1:num2
    i1 = i;
    i2 = mod(i, num2) + 1;
    edge2 = diff(vtx2([i1, i2], :));

    % Find out which point of vtx1 lies inner side,
    % and deal with edge crossing.
    k = 1;
    for j = 1:num1
        j1 = j;
        j2 = mod(j, num1) + 1;
        z1 = det([edge2; vtx1(j1, :) - vtx2(i2, :)]);
        z2 = det([edge2; vtx1(j2, :) - vtx2(i2, :)]);

        if z1 * z0 > -eps_z && z2 * z0 > -eps_z
            % Two inner (at edge) points. Record the first one
            vtx(k, :) = vtx1(j1, :);
            k = k + 1;
        elseif z1 * z0 > -eps_z && z2 * z0 < eps_z
            % Inner (at edge) --> outer. Record the first point and the intersection point.
            vtx(k, :) = vtx1(j1, :);
            s = (vtx2(i1, :) - vtx1(j1, :)) / [diff(vtx2([i1, i2], :)); diff(vtx1([j1, j2], :))];
            vtx(k + 1, :) = s(2) * diff(vtx1([j1, j2], :)) + vtx1(j1, :);
            k = k + 2;
        elseif z1 * z0 < eps_z && z2 * z0 > -eps_z
            % Outer --> inner (at edge). Record the intersection point.
            s = (vtx2(i1, :) - vtx1(j1, :)) / [diff(vtx2([i1, i2], :)); diff(vtx1([j1, j2], :))];
            vtx(k, :) = s(2) * diff(vtx1([j1, j2], :)) + vtx1(j1, :);
            k = k + 1;
        end
    end
    num1 = k - 1;
    vtx1 = vtx(1:num1, :);
end
vtx = vtx1;
end
