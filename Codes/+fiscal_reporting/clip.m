function values = clip(values, limit)
% Preserve deliberate infinity reports while excluding large finite ratios.
values(isfinite(values) & abs(values) > limit) = NaN;
end
