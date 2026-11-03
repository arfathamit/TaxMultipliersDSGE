function value = percentile(sample, probabilities)
% MATLAB prctile's midpoint interpolation, retaining infinite endpoints.
% Avoid Inf-Inf arithmetic, which otherwise turns all-Inf intervals into NaN.
sample = sort(sample(~isnan(sample)));
value = NaN(size(probabilities));
if isempty(sample), return; end
if any(probabilities(:) < 0 | probabilities(:) > 100)
    error('fiscal_reporting:InvalidPercentile', 'Percentiles must be in [0,100].');
end
n = numel(sample);
for k = 1:numel(probabilities)
    rank = min(n, max(1, n * probabilities(k)/100 + 0.5));
    lo = floor(rank); hi = ceil(rank);
    a = sample(lo); b = sample(hi);
    if lo == hi || a == b
        value(k) = a;
    elseif isinf(a) && isinf(b)
        value(k) = NaN;  % Interpolation between -Inf and +Inf is undefined.
    elseif isinf(a)
        value(k) = a;
    elseif isinf(b)
        value(k) = b;
    else
        value(k) = a + (rank-lo)*(b-a);
    end
end
end
