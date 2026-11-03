function value = multiplier(output, own_cost, net_cost)
% Reporting convention for expansionary fiscal experiments.
% OUTPUT is PV GDP change in the tax-cut/incentive-increase direction.
% OWN_COST is the signed denominator used by the requested normalization.
% NET_COST is minus PV total tax revenue change in that same direction.
% Pass NaN for NET_COST for spending/transfer experiments (no override).
% Positive output with zero or negative net revenue cost is reported +Inf.
% This is a reporting convention, not the algebraic value of OUTPUT/OWN_COST.
if ~isequal(size(output), size(own_cost), size(net_cost))
    error('fiscal_reporting:SizeMismatch', 'Inputs must have identical sizes.');
end
tol = 1e-10;
value = NaN(size(output));
valid = isfinite(output) & isfinite(own_cost) & abs(own_cost) > tol;
value(valid) = output(valid) ./ own_cost(valid);
self_financing = isfinite(output) & isfinite(net_cost) ...
    & output > tol & net_cost <= tol;
value(self_financing) = Inf;
end
