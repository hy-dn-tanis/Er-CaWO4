% normalize any quantity to the mean along frequency and field
function [Z_norm] = normalize_freq_field_mean(Z)

f_mean = zeros([1, size(Z, 2)]);
B_mean = zeros([1, size(Z, 1)]);
Z_sub_f = zeros(size(Z));
Z_sub_f_sub_B = zeros(size(Z));

for bb = 1:size(Z, 2) % subtract mean along frequency
    f_mean(bb) = mean(Z(:, bb));
    Z_sub_f(:, bb) = Z(:, bb) - f_mean(bb);
end

for ff = 1:size(Z, 1) % subtract mean along field
    B_mean(ff) = mean(Z_sub_f(ff, :));
    Z_sub_f_sub_B(ff, :) = Z_sub_f(ff, :) - B_mean(ff);
end

Z_norm = Z_sub_f_sub_B;

end