function axis_pdf = make_axis_pdf(zenith_dist, roll_dist)
% Make axis distribution PDF.
%
% INPUT
%   zenith_dist:            [mean, sigma] or []
%   roll_dist:              [mean, sigma] or []

if nargin == 0 || isempty(zenith_dist)
    % Uniform distribution
    a_pdf = @(llr) ones(size(llr, 1), 1) / (180 * 180 * 4 / pi);
else
    % Gaussian distribution
    zenith = zenith_dist(1:2);
    tmp_x = linspace(-90, 90, 50000);
    tmp_pdf = exp(- (90 - tmp_x - zenith(1)).^2/2 / zenith(2)^2) / zenith(2);
    zen_total = (sum(tmp_pdf) * (tmp_x(2) - tmp_x(1)));
    a_pdf = @(llr) (exp(- (90 - asind(sind(llr(:, 2))) - zenith(1)).^2/2 / zenith(2)^2) / ...
        zenith(2) / zen_total);
end

if nargin == 1 || isempty(roll_dist)
    r_pdf = @(llr) ones(size(llr, 1), 1) / 360;
else
    roll = roll_dist(1:2);
    tmp_x = linspace(0, 360, 50000);
    tmp_pdf = exp(- (tmp_x - roll(1)).^2/2 / roll(2)^2) / roll(2);
    zen_total = (sum(tmp_pdf) * (tmp_x(2) - tmp_x(1)));
    r_pdf = @(llr) (exp(- (llr(:, 3) - roll(1)).^2/2 / roll(2)^2) / ...
        roll(2) / zen_total);
end

axis_pdf = @(llr) a_pdf(llr) .* r_pdf(llr);
end
