function res = check_coplanar(pts, tol)
% Check if all points are coplanar
%
% INPUT
%   pts:    n*d, n points
%   tol:    scalar, tolerance for rank
%
% OUTPUT
%   res:    true or false

num = size(pts, 1);
if nargin == 1
    tol = 1e-10;
end

if num <= 3
    res = true;
    return;
end

e1 = pts(2, :) - pts(1, :);
e2 = e1;
i = 3;
while rank([e1; e2], tol) < 2 && i <= num
    e2 = pts(i, :) - pts(1, :);
    i = i + 1;
end
if i > num
    res = true;
    return;
end

d = bsxfun(@minus, pts, pts(1, :));
res = true;
while i <= num
    if rank([e1; e2; d(i, :)], tol) > 2
        res = false;
        return;
    end
    i = i + 1;
end
end
