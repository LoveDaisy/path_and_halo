function axis_pdf = make_axis_pdf(zenith, roll)
% Make axis distribution PDF.
%
% INPUT
%   zenith:         [mean, sigma] or []
%   roll:           [mean, sigma] or []
%
% OUTPUT
%   axis_pdf:       function handle. It has the form: p = fun(llr).
%                   It can be decomposed into two parts, [longitude, latitude] is direction part,
%                   [roll] is rotation part. These two part can be integrated separately.
%                       integrate(pdf_ll, d_lat, d_lon) * integrate(pdf_r, d_roll) = 1
%                   So we can normalize them to 1 respectively, and take the total integration constant
%                   at last.

% First part, axis direction.
% This PDF represents the density of axis direction, i.e. the probability per unit area on sphere.
% It is normalized so that C * integrate(pdf_ll, d_area) = 1.
if nargin == 0 || isempty(zenith)
    % Uniform distribution
    a_pdf = @(llr) ones(size(llr, 1), 1) / (4 * pi);
else
    % Gaussian distribution
    tmp_x = linspace(-90, 90, 50000);
    tmp_pdf = exp(- (90 - tmp_x - zenith(1)).^2/2 / zenith(2)^2) / zenith(2);
    zen_total = (sum(tmp_pdf .* cosd(tmp_x)) * (tmp_x(2) - tmp_x(1)) * pi / 180 * 2 * pi);
    a_pdf = @(llr) (exp(- (90 - asind(sind(llr(:, 2))) - zenith(1)).^2/2 / zenith(2)^2) / ...
        zenith(2) / zen_total);
end

% Second part, roll.
% This PDF represents the density of rolling, i.e. the probability per unit rolling angle.
% It is normalized so that C * integrate(pdf_roll, d_theta) = 1.
if nargin <= 1 || isempty(roll)
    r_pdf = @(llr) ones(size(llr, 1), 1) / (2 * pi);
else
    tmp_x = linspace(0, 360, 50000);
    tmp_pdf = exp(- (tmp_x - roll(1)).^2/2 / roll(2)^2) / roll(2);
    zen_total = (sum(tmp_pdf) * (tmp_x(2) - tmp_x(1)) * pi / 180);
    r_pdf = @(llr) (exp(- (llr(:, 3) - roll(1)).^2/2 / roll(2)^2) / ...
        roll(2) / zen_total);
end

axis_pdf = @(llr) a_pdf(llr) .* r_pdf(llr);
end
