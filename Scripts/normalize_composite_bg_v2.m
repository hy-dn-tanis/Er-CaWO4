% normalize any quantity to a composite background
function [Z_norm, bg_fields] = normalize_composite_bg_v2(Z, f, B, bg_limits)
% Z: raw data to normalize; dimensions (frequency points)x(field points)
% f: frequency points (GHz)
% B: field points (mT)
% bg_limits: [frequencies (GHz), fields (mT)] starting points of intervals for the spectra patches
% Z_norm: normalized data; dimensions (frequency points)x(field points)
% bg_fields: array of fields of the patches; for plotting the sections from which the background was taken from

bg = [];
bg_fields = zeros(size(f));

ind2 = 0;
% concatenate spectra for each interval at different fields
for i = 1:size(bg_limits, 1) - 1
    ind1 = ind2 + 1; % start index of interval
    ind2 = find(f/1e9 < bg_limits(i+1, 1), 1, 'last'); % stop index of interval
    ind_B = find(B >= bg_limits(i, 2), 1, 'first'); % field index of interval
    bg = [bg transpose(Z(ind1:ind2, ind_B))];
    bg_fields(ind1:ind2) = bg_limits(i, 2);
end

% append the last point from the same field as the last interval
bg = transpose([bg Z(end, ind_B)]);
bg_fields(end) = bg_limits(end-1, 2);
Z_norm = Z - bg;

end



