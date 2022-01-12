function axis_pdf = make_axis_pdf(zenith_dist, roll_dist)
% Make axis distribution PDF.
%
% INPUT
%   zenith_dist:            [type, mean, sigma]
%   roll_dist:              [type, mean, sigma]

p = inputParser;
p.addRequired('zenith_dist', @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 3}));
p.addRequired('roll_dist', @(x) isempty(x) || validateattributes(x, {'numeric'}, {'vector', 'numel', 3}));
p.parse(zenith_dist, roll_dist);

if zenith_dist(1) == 0
    % Uniform distribution
    a_pdf = @(llr) ones(size(llr, 1), 1) / (180 * 180 * 4 / pi);
elseif zenith_dist(1) == 1
    % Gaussian distribution
    zenith = zenith_dist(2:3);
    tmp_x = linspace(-90, 90, 50000);
    tmp_pdf = exp(- (90 - tmp_x - zenith(1)).^2/2 / zenith(2)^2) / zenith(2);
    zen_total = (sum(tmp_pdf) * (tmp_x(2) - tmp_x(1)));
    a_pdf = @(llr) (exp(- (90 - asind(sind(llr(:, 2))) - zenith(1)).^2/2 / zenith(2)^2) / ...
        zenith(2) / zen_total);
end

if isempty(roll_dist) || roll_dist(1) == 0
    r_pdf = @(llr) ones(size(llr, 1), 1) / 360;
elseif roll_dist(1) == 1
    roll = roll_dist(2:3);
    tmp_x = linspace(0, 360, 50000);
    tmp_pdf = exp(- (tmp_x - roll(1)).^2/2 / roll(2)^2) / roll(2);
    zen_total = (sum(tmp_pdf) * (tmp_x(2) - tmp_x(1)));
    r_pdf = @(llr) (exp(- (llr(:, 3) - roll(1)).^2/2 / roll(2)^2) / ...
        roll(2) / zen_total);
end

axis_pdf = @(llr) a_pdf(llr) .* r_pdf(llr);
end
