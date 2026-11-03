function test_fiscal_reporting()
% Run in MATLAB or Octave: addpath('Codes/tests'); test_fiscal_reporting
addpath(fileparts(fileparts(mfilename('fullpath'))));

% A positive GDP effect is self-financing based on TOTAL revenue, even
% when the instrument's own outlay (the usual denominator) stays positive.
output = [2 2 2 -2 0 NaN 2];
own = [4 4 4 4 4 4 4];
net = [1 -1 0 -1 -1 -1 NaN];
r = fiscal_reporting.multiplier(output, own, net);
assert(r(1)==0.5 && all(r(2:3)==Inf));
assert(r(4)==-0.5 && r(5)==0 && isnan(r(6)) && r(7)==0.5);
assert(isnan(fiscal_reporting.multiplier(0,0,0)));

% Positive innovations have opposite meanings for tax rates and credits.
% A rate increase lowers Y and lowers revenue: reversing it is a tax cut
% that raises BOTH Y and revenue, hence +Inf under the reporting convention.
raw_y = -2; raw_rev = -1; expansion_direction = -1;
assert(fiscal_reporting.multiplier(expansion_direction*raw_y,raw_rev, ...
    -expansion_direction*raw_rev)==Inf);
% Credit increase raises both Y and revenue; the same rule applies.
assert(fiscal_reporting.multiplier(2,4,-1)==Inf);
% A contractionary experiment is not relabeled as a beneficial free policy.
assert(fiscal_reporting.multiplier(-2,4,-1)==-0.5);

% Classification, clipping, and summaries all retain the infinite reports.
clipped = fiscal_reporting.clip([2 51 -51 Inf NaN],50);
assert(clipped(1)==2 && all(isnan(clipped(2:3))) ...
    && clipped(4)==Inf && isnan(clipped(5)));
assert(all(fiscal_reporting.percentile([Inf Inf NaN],[5 50 95])==Inf));
assert(isequal(fiscal_reporting.percentile([1 2 3 Inf Inf],[5 50 95]),[1 3 Inf]));
assert(fiscal_reporting.percentile([1 Inf],50)==Inf);
assert(all(isnan(fiscal_reporting.percentile([NaN NaN],[5 50 95]))));
finite = [2 -3 0 7 4 1 9 12 1 5];
p = [0 5 25 50 75 95 100];
assert(max(abs(fiscal_reporting.percentile(finite,p)-prctile(finite,p)))<1e-12);

% Existing saved-mode examples: revenue gain makes the ITC +Inf, while
% the negative-output expensing result retains its negative finite value.
assert(fiscal_reporting.multiplier(2.259602913,.922312833,-.132980791)==Inf);
e = fiscal_reporting.multiplier(-.200046321,.622326747,.227497015);
assert(abs(e-(-.321449017))<1e-8);
% Multidimensional storage used by post_estimation.
v = reshape([1 2 -1 3], [1 2 1 2]);
cost = ones(size(v)); net = reshape([-1 1 -1 0],size(v));
r = fiscal_reporting.multiplier(v,cost,net);
assert(isequal(size(r),size(v)) && isequal(r(:)',[Inf 2 -1 Inf]));
fprintf('Fiscal reporting tests passed.\n');
end
