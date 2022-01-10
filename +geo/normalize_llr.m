function llr = normalize_llr(llr)
% Normalize LLR representation.
% Longitude is in range [-180, 180], latitude is in range [-90, 90], and roll is
% in range [0, 360]

llr(:, 1) = mod(llr(:, 1), 360);
idx = llr(:, 1) > 180;
llr(idx, 1) = llr(idx, 1) - 360;

llr(:, 2) = mod(llr(:, 2), 180);
idx = llr(:, 2) > 90;
llr(idx, 2) = llr(idx, 2) - 180;

llr(:, 3) = mod(llr(:, 3), 360);
end