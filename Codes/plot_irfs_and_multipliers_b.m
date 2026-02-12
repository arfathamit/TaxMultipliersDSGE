%% plot_irfs_and_multipliers_b.m
%% =========================================================================
%% Generate IRF plots and multiplier tables
%% Replicates presentation style of Mertens & Ravn (2013) and Zubairy (2014)
%% 
%% VERSION B: For model with 6 TAX SHOCKS (no government spending shock)
%%   - tau_c (consumption tax)
%%   - tau_l (labor income tax)
%%   - tau_k (capital income tax)
%%   - tau_i (interest income tax)
%%   - tau_d (dividend tax)
%%   - tau_ic (investment tax credit)
%%
%% Run this AFTER running Dynare on tax_dsge_fiscal_shocks_b.mod
%% =========================================================================

close all;

%% =========================================================================
%% PART 1: SETUP AND DATA EXTRACTION
%% =========================================================================

% Check if Dynare results exist
if ~exist('oo_', 'var')
    error('Dynare results not found. Run dynare tax_dsge_fiscal_shocks_b first.');
end

fprintf('\n');
fprintf('================================================================\n');
fprintf('  FISCAL MULTIPLIER ANALYSIS - Post-Processing Script (Option B)\n');
fprintf('================================================================\n');
fprintf('\n');

% Define horizon
T = 40;  % quarters
quarters = 0:T-1;

% Horizons for multiplier table (matching Zubairy 2014)
horizons = [1, 4, 12, 20];
n_horizons = length(horizons);

%% =========================================================================
%% PART 2: EXTRACT STEADY STATE VALUES
%% =========================================================================

fprintf('Extracting steady-state values...\n');

% Function to get steady state value by variable name
get_ss = @(varname) oo_.steady_state(strmatch(varname, M_.endo_names, 'exact'));

% Get steady state values
y_ss = get_ss('y');
c_ss = get_ss('c');
i_ss = get_ss('i');
g_ss = get_ss('g');
l_ss = get_ss('l');
w_ss = get_ss('w');
r_k_ss = get_ss('r_k');
k_ss = get_ss('k');
d_ss = get_ss('d');
b_ss = get_ss('b');
r_ss = get_ss('r');

% Additional steady state values needed for tax revenue computations
p_ss = get_ss('p');
v_ss = get_ss('v');
k_tau_ss = get_ss('k_tau');

% Get steady-state tax rates
tau_c_ss = get_ss('tau_c');
tau_l_ss = get_ss('tau_l');
tau_k_ss = get_ss('tau_k');
tau_i_ss = get_ss('tau_i');
tau_d_ss = get_ss('tau_d');
tau_ic_ss = get_ss('tau_ic');

% Get discount factor
beta = M_.params(strmatch('betta', M_.param_names, 'exact'));
delta_tau = M_.params(strmatch('delta_tau', M_.param_names, 'exact'));
e_tau     = M_.params(strmatch('e_tau', M_.param_names, 'exact'));

R = 1/beta - 1;  % Discount rate for present value calculations

fprintf('  GDP (y_ss) = %.4f\n', y_ss);
fprintf('  Consumption (c_ss) = %.4f\n', c_ss);
fprintf('  Investment (i_ss) = %.4f\n', i_ss);
fprintf('  Gov spending (g_ss) = %.4f\n', g_ss);
fprintf('  Discount rate (R) = %.4f\n', R);
fprintf('\n');

%% =========================================================================
%% PART 3: CALCULATE TAX REVENUE BASES (for normalization)
%% =========================================================================

fprintf('Calculating tax revenue bases for normalization...\n');

% Tax revenue bases (steady state)
% These are used to convert tax rate changes to revenue changes
T_c_base  = c_ss;                                               % Consumption tax base
T_l_base  = (w_ss * l_ss) / p_ss;                                 % Labor income tax base (real wage bill)
T_k_base  = (r_k_ss * v_ss * k_ss) / p_ss - (delta_tau * k_tau_ss) - (e_tau * i_ss); % Capital income tax base (net of allowances)
T_i_base  = (r_ss * b_ss) / p_ss;                                 % Interest income tax base (real interest income)
T_d_base  = d_ss / p_ss;                                          % Dividend tax base (real dividends)
T_ic_base = i_ss;                                                 % ITC base (investment)

% Steady-state fiscal flows
Rev_c_ss   = tau_c_ss  * T_c_base;
Rev_l_ss   = tau_l_ss  * T_l_base;
Rev_k_ss   = tau_k_ss  * T_k_base;
Rev_i_ss   = tau_i_ss  * T_i_base;
Rev_d_ss   = tau_d_ss  * T_d_base;
Spend_ic_ss = tau_ic_ss * T_ic_base;                               % ITC subsidy outlay
fprintf('  Consumption tax revenue = %.4f (%.1f%% of GDP)\n', Rev_c_ss, 100*Rev_c_ss/y_ss);
fprintf('  Labor tax revenue = %.4f (%.1f%% of GDP)\n', Rev_l_ss, 100*Rev_l_ss/y_ss);
fprintf('  Capital tax revenue = %.4f (%.1f%% of GDP)\n', Rev_k_ss, 100*Rev_k_ss/y_ss);
fprintf('  Interest tax revenue = %.4f (%.1f%% of GDP)\n', Rev_i_ss, 100*Rev_i_ss/y_ss);
fprintf('  Dividend tax revenue = %.4f (%.1f%% of GDP)\n', Rev_d_ss, 100*Rev_d_ss/y_ss);
fprintf('  ITC subsidy outlay = %.4f (%.1f%% of GDP)\n', Spend_ic_ss, 100*Spend_ic_ss/y_ss);
fprintf('\n');

%% =========================================================================
%% PART 4: EXTRACT ALL IRFs
%% =========================================================================

fprintf('Extracting impulse response functions...\n');

% --- Consumption Tax Shock ---
irf_y_tc = oo_.irfs.y_e_tau_c;
irf_c_tc = oo_.irfs.c_e_tau_c;
irf_i_tc = oo_.irfs.i_e_tau_c;
irf_l_tc = oo_.irfs.l_e_tau_c;
irf_tc_tc = oo_.irfs.tau_c_e_tau_c;
irf_p_tc = oo_.irfs.p_e_tau_c;
irf_w_tc = oo_.irfs.w_e_tau_c;
irf_r_tc = oo_.irfs.r_e_tau_c;
irf_rk_tc = oo_.irfs.r_k_e_tau_c;
irf_v_tc = oo_.irfs.v_e_tau_c;
irf_k_tc = oo_.irfs.k_e_tau_c;
irf_ktau_tc = oo_.irfs.k_tau_e_tau_c;
irf_b_tc = oo_.irfs.b_e_tau_c;
irf_d_tc = oo_.irfs.d_e_tau_c;
irf_rev_c_tc = oo_.irfs.rev_c_e_tau_c;
irf_rev_c_tc = oo_.irfs.rev_c_e_tau_c;
irf_rev_l_tc = oo_.irfs.rev_l_e_tau_c;
irf_rev_k_tc = oo_.irfs.rev_k_e_tau_c;
irf_rev_i_tc = oo_.irfs.rev_i_e_tau_c;
irf_rev_d_tc = oo_.irfs.rev_d_e_tau_c;
irf_rev_total_tc = irf_rev_c_tc + irf_rev_l_tc + irf_rev_k_tc + irf_rev_i_tc + irf_rev_d_tc;



% --- Labor Tax Shock ---
irf_y_tl = oo_.irfs.y_e_tau_l;
irf_c_tl = oo_.irfs.c_e_tau_l;
irf_i_tl = oo_.irfs.i_e_tau_l;
irf_l_tl = oo_.irfs.l_e_tau_l;
irf_tl_tl = oo_.irfs.tau_l_e_tau_l;
irf_p_tl = oo_.irfs.p_e_tau_l;
irf_w_tl = oo_.irfs.w_e_tau_l;
irf_r_tl = oo_.irfs.r_e_tau_l;
irf_rk_tl = oo_.irfs.r_k_e_tau_l;
irf_v_tl = oo_.irfs.v_e_tau_l;
irf_k_tl = oo_.irfs.k_e_tau_l;
irf_ktau_tl = oo_.irfs.k_tau_e_tau_l;
irf_b_tl = oo_.irfs.b_e_tau_l;
irf_d_tl = oo_.irfs.d_e_tau_l;
irf_rev_l_tl = oo_.irfs.rev_l_e_tau_l;
irf_rev_c_tl = oo_.irfs.rev_c_e_tau_l;
irf_rev_l_tl = oo_.irfs.rev_l_e_tau_l;
irf_rev_k_tl = oo_.irfs.rev_k_e_tau_l;
irf_rev_i_tl = oo_.irfs.rev_i_e_tau_l;
irf_rev_d_tl = oo_.irfs.rev_d_e_tau_l;
irf_rev_total_tl = irf_rev_c_tl + irf_rev_l_tl + irf_rev_k_tl + irf_rev_i_tl + irf_rev_d_tl;



% --- Capital Tax Shock ---
irf_y_tk = oo_.irfs.y_e_tau_k;
irf_c_tk = oo_.irfs.c_e_tau_k;
irf_i_tk = oo_.irfs.i_e_tau_k;
irf_l_tk = oo_.irfs.l_e_tau_k;
irf_tk_tk = oo_.irfs.tau_k_e_tau_k;
irf_p_tk = oo_.irfs.p_e_tau_k;
irf_w_tk = oo_.irfs.w_e_tau_k;
irf_r_tk = oo_.irfs.r_e_tau_k;
irf_rk_tk = oo_.irfs.r_k_e_tau_k;
irf_v_tk = oo_.irfs.v_e_tau_k;
irf_k_tk = oo_.irfs.k_e_tau_k;
irf_ktau_tk = oo_.irfs.k_tau_e_tau_k;
irf_b_tk = oo_.irfs.b_e_tau_k;
irf_d_tk = oo_.irfs.d_e_tau_k;
irf_rev_k_tk = oo_.irfs.rev_k_e_tau_k;
irf_rev_c_tk = oo_.irfs.rev_c_e_tau_k;
irf_rev_l_tk = oo_.irfs.rev_l_e_tau_k;
irf_rev_k_tk = oo_.irfs.rev_k_e_tau_k;
irf_rev_i_tk = oo_.irfs.rev_i_e_tau_k;
irf_rev_d_tk = oo_.irfs.rev_d_e_tau_k;
irf_rev_total_tk = irf_rev_c_tk + irf_rev_l_tk + irf_rev_k_tk + irf_rev_i_tk + irf_rev_d_tk;



% --- Interest Tax Shock ---
irf_y_ti = oo_.irfs.y_e_tau_i;
irf_c_ti = oo_.irfs.c_e_tau_i;
irf_i_ti = oo_.irfs.i_e_tau_i;
irf_l_ti = oo_.irfs.l_e_tau_i;
irf_ti_ti = oo_.irfs.tau_i_e_tau_i;
irf_p_ti = oo_.irfs.p_e_tau_i;
irf_w_ti = oo_.irfs.w_e_tau_i;
irf_r_ti = oo_.irfs.r_e_tau_i;
irf_rk_ti = oo_.irfs.r_k_e_tau_i;
irf_v_ti = oo_.irfs.v_e_tau_i;
irf_k_ti = oo_.irfs.k_e_tau_i;
irf_ktau_ti = oo_.irfs.k_tau_e_tau_i;
irf_b_ti = oo_.irfs.b_e_tau_i;
irf_d_ti = oo_.irfs.d_e_tau_i;
irf_rev_i_ti = oo_.irfs.rev_i_e_tau_i;
irf_rev_c_ti = oo_.irfs.rev_c_e_tau_i;
irf_rev_l_ti = oo_.irfs.rev_l_e_tau_i;
irf_rev_k_ti = oo_.irfs.rev_k_e_tau_i;
irf_rev_i_ti = oo_.irfs.rev_i_e_tau_i;
irf_rev_d_ti = oo_.irfs.rev_d_e_tau_i;
irf_rev_total_ti = irf_rev_c_ti + irf_rev_l_ti + irf_rev_k_ti + irf_rev_i_ti + irf_rev_d_ti;



% --- Dividend Tax Shock ---
irf_y_td = oo_.irfs.y_e_tau_d;
irf_c_td = oo_.irfs.c_e_tau_d;
irf_i_td = oo_.irfs.i_e_tau_d;
irf_l_td = oo_.irfs.l_e_tau_d;
irf_td_td = oo_.irfs.tau_d_e_tau_d;
irf_p_td = oo_.irfs.p_e_tau_d;
irf_w_td = oo_.irfs.w_e_tau_d;
irf_r_td = oo_.irfs.r_e_tau_d;
irf_rk_td = oo_.irfs.r_k_e_tau_d;
irf_v_td = oo_.irfs.v_e_tau_d;
irf_k_td = oo_.irfs.k_e_tau_d;
irf_ktau_td = oo_.irfs.k_tau_e_tau_d;
irf_b_td = oo_.irfs.b_e_tau_d;
irf_d_td = oo_.irfs.d_e_tau_d;
irf_rev_d_td = oo_.irfs.rev_d_e_tau_d;
irf_rev_c_td = oo_.irfs.rev_c_e_tau_d;
irf_rev_l_td = oo_.irfs.rev_l_e_tau_d;
irf_rev_k_td = oo_.irfs.rev_k_e_tau_d;
irf_rev_i_td = oo_.irfs.rev_i_e_tau_d;
irf_rev_d_td = oo_.irfs.rev_d_e_tau_d;
irf_rev_total_td = irf_rev_c_td + irf_rev_l_td + irf_rev_k_td + irf_rev_i_td + irf_rev_d_td;



% --- Investment Tax Credit Shock ---
irf_y_tic = oo_.irfs.y_e_tau_ic;
irf_c_tic = oo_.irfs.c_e_tau_ic;
irf_i_tic = oo_.irfs.i_e_tau_ic;
irf_l_tic = oo_.irfs.l_e_tau_ic;
irf_tic_tic = oo_.irfs.tau_ic_e_tau_ic;
irf_p_tic = oo_.irfs.p_e_tau_ic;
irf_w_tic = oo_.irfs.w_e_tau_ic;
irf_r_tic = oo_.irfs.r_e_tau_ic;
irf_rk_tic = oo_.irfs.r_k_e_tau_ic;
irf_v_tic = oo_.irfs.v_e_tau_ic;
irf_k_tic = oo_.irfs.k_e_tau_ic;
irf_ktau_tic = oo_.irfs.k_tau_e_tau_ic;
irf_b_tic = oo_.irfs.b_e_tau_ic;
irf_d_tic = oo_.irfs.d_e_tau_ic;
irf_spend_ic_tic = oo_.irfs.spend_ic_e_tau_ic;
irf_rev_c_tic = oo_.irfs.rev_c_e_tau_ic;
irf_rev_l_tic = oo_.irfs.rev_l_e_tau_ic;
irf_rev_k_tic = oo_.irfs.rev_k_e_tau_ic;
irf_rev_i_tic = oo_.irfs.rev_i_e_tau_ic;
irf_rev_d_tic = oo_.irfs.rev_d_e_tau_ic;
irf_rev_total_tic = irf_rev_c_tic + irf_rev_l_tic + irf_rev_k_tic + irf_rev_i_tic + irf_rev_d_tic;
irf_spend_ic_tic = oo_.irfs.spend_ic_e_tau_ic;



fprintf('  Extracted IRFs for all 6 tax shocks.\n');
%% =========================================================================
%% PART 4B: NORMALIZE SHOCKS (Zubairy-style)
%% =========================================================================
% Tax shocks are normalized so that the implied impact change in total tax revenue
% equals -1% of steady-state total tax revenue (tax cut convention).
% The ITC shock is normalized so that the implied impact change in ITC outlays
% equals +1% of steady-state GDP (since steady-state ITC outlays can be zero).

Rev_total_ss = Rev_c_ss + Rev_l_ss + Rev_k_ss + Rev_i_ss + Rev_d_ss;
target_tax_rev_drop = 0.01 * Rev_total_ss;
target_itc_outlay   = 0.01 * y_ss;

% Guard against zero impact responses
if abs(irf_rev_total_tc(1)) < 1e-12, error('Impact total revenue response is zero for tau_c shock.'); end
if abs(irf_rev_total_tl(1)) < 1e-12, error('Impact total revenue response is zero for tau_l shock.'); end
if abs(irf_rev_total_tk(1)) < 1e-12, error('Impact total revenue response is zero for tau_k shock.'); end
if abs(irf_rev_total_ti(1)) < 1e-12, error('Impact total revenue response is zero for tau_i shock.'); end
if abs(irf_rev_total_td(1)) < 1e-12, error('Impact total revenue response is zero for tau_d shock.'); end
if abs(irf_spend_ic_tic(1)) < 1e-12, error('Impact ITC outlay response is zero for tau_ic shock.'); end

% Scaling factors for tax cuts: apply a negative sign so the normalized impulse is a cut
scale_tc = target_tax_rev_drop / irf_rev_total_tc(1);
scale_tl = target_tax_rev_drop / irf_rev_total_tl(1);
scale_tk = target_tax_rev_drop / irf_rev_total_tk(1);
scale_ti = target_tax_rev_drop / irf_rev_total_ti(1);
scale_td = target_tax_rev_drop / irf_rev_total_td(1);

% Scaling factor for ITC increase (outlay-based normalization)
scale_tic = target_itc_outlay / irf_spend_ic_tic(1);

% Normalized IRFs (already in the reporting direction: tax cuts and ITC increase)
irf_y_tc_n = -scale_tc * irf_y_tc;
irf_y_tl_n = -scale_tl * irf_y_tl;
irf_y_tk_n = -scale_tk * irf_y_tk;
irf_y_ti_n = -scale_ti * irf_y_ti;
irf_y_td_n = -scale_td * irf_y_td;
irf_y_tic_n =  scale_tic * irf_y_tic;

irf_c_tc_n = -scale_tc * irf_c_tc;
irf_c_tl_n = -scale_tl * irf_c_tl;
irf_c_tk_n = -scale_tk * irf_c_tk;
irf_c_ti_n = -scale_ti * irf_c_ti;
irf_c_td_n = -scale_td * irf_c_td;
irf_c_tic_n =  scale_tic * irf_c_tic;

irf_i_tc_n = -scale_tc * irf_i_tc;
irf_i_tl_n = -scale_tl * irf_i_tl;
irf_i_tk_n = -scale_tk * irf_i_tk;
irf_i_ti_n = -scale_ti * irf_i_ti;
irf_i_td_n = -scale_td * irf_i_td;
irf_i_tic_n =  scale_tic * irf_i_tic;

irf_l_tc_n = -scale_tc * irf_l_tc;
irf_l_tl_n = -scale_tl * irf_l_tl;
irf_l_tk_n = -scale_tk * irf_l_tk;
irf_l_ti_n = -scale_ti * irf_l_ti;
irf_l_td_n = -scale_td * irf_l_td;
irf_l_tic_n =  scale_tic * irf_l_tic;

irf_rev_total_tc_n = -scale_tc * irf_rev_total_tc;
irf_rev_total_tl_n = -scale_tl * irf_rev_total_tl;
irf_rev_total_tk_n = -scale_tk * irf_rev_total_tk;
irf_rev_total_ti_n = -scale_ti * irf_rev_total_ti;
irf_rev_total_td_n = -scale_td * irf_rev_total_td;

irf_spend_ic_tic_n = scale_tic * irf_spend_ic_tic;

% Also normalize tax rate paths (useful for plots)
irf_tc_tc_n = -scale_tc * irf_tc_tc;
irf_tl_tl_n = -scale_tl * irf_tl_tl;
irf_tk_tk_n = -scale_tk * irf_tk_tk;
irf_ti_ti_n = -scale_ti * irf_ti_ti;
irf_td_td_n = -scale_td * irf_td_td;
irf_tic_tic_n = scale_tic * irf_tic_tic;

fprintf('Normalization complete: tax cuts imply -1%% impact total tax revenue, ITC implies +1%% GDP impact outlay.\n\n');

fprintf('\n');

%% =========================================================================
%% PART 4C: RAW IRFS FOR PLOTS (NO NORMALIZATION)
%% =========================================================================
% For plots, keep the original Dynare shock size.
% Use the same sign convention as the multiplier reporting:
%   - Tax shocks are shown as TAX CUTS (flip sign)
%   - ITC shock is shown as an INCREASE in the credit (no sign flip)

%% =========================================================================
%% Plot scaling option
%% If true, IRF plots use the normalized IRFs (same scaling used for multipliers).
%% If false, IRF plots use raw Dynare IRFs (with sign flips for tax cuts only).
SCALE_PLOTS_TO_NORMALIZATION = true;

irf_y_tc_p   = -irf_y_tc;
irf_y_tl_p   = -irf_y_tl;
irf_y_tk_p   = -irf_y_tk;
irf_y_ti_p   = -irf_y_ti;
irf_y_td_p   = -irf_y_td;
irf_y_tic_p  =  irf_y_tic;

irf_c_tc_p   = -irf_c_tc;
irf_c_tl_p   = -irf_c_tl;
irf_c_tk_p   = -irf_c_tk;
irf_c_ti_p   = -irf_c_ti;
irf_c_td_p   = -irf_c_td;
irf_c_tic_p  =  irf_c_tic;

irf_i_tc_p   = -irf_i_tc;
irf_i_tl_p   = -irf_i_tl;
irf_i_tk_p   = -irf_i_tk;
irf_i_ti_p   = -irf_i_ti;
irf_i_td_p   = -irf_i_td;
irf_i_tic_p  =  irf_i_tic;

irf_l_tc_p   = -irf_l_tc;
irf_l_tl_p   = -irf_l_tl;
irf_l_tk_p   = -irf_l_tk;
irf_l_ti_p   = -irf_l_ti;
irf_l_td_p   = -irf_l_td;
irf_l_tic_p  =  irf_l_tic;

% Tax rate paths (percentage points), plotted as cuts for taxes
irf_tc_tc_p  = -irf_tc_tc;
irf_tl_tl_p  = -irf_tl_tl;
irf_tk_tk_p  = -irf_tk_tk;
irf_ti_ti_p  = -irf_ti_ti;
irf_td_td_p  = -irf_td_td;
irf_tic_tic_p = irf_tic_tic;

if SCALE_PLOTS_TO_NORMALIZATION
    irf_c_tc_p = irf_c_tc_n;
    irf_c_td_p = irf_c_td_n;
    irf_c_ti_p = irf_c_ti_n;
    irf_c_tic_p = irf_c_tic_n;
    irf_c_tk_p = irf_c_tk_n;
    irf_c_tl_p = irf_c_tl_n;
    irf_i_tc_p = irf_i_tc_n;
    irf_i_td_p = irf_i_td_n;
    irf_i_ti_p = irf_i_ti_n;
    irf_i_tic_p = irf_i_tic_n;
    irf_i_tk_p = irf_i_tk_n;
    irf_i_tl_p = irf_i_tl_n;
    irf_l_tc_p = irf_l_tc_n;
    irf_l_td_p = irf_l_td_n;
    irf_l_ti_p = irf_l_ti_n;
    irf_l_tic_p = irf_l_tic_n;
    irf_l_tk_p = irf_l_tk_n;
    irf_l_tl_p = irf_l_tl_n;
    irf_tc_tc_p = irf_tc_tc_n;
    irf_td_td_p = irf_td_td_n;
    irf_ti_ti_p = irf_ti_ti_n;
    irf_tic_tic_p = irf_tic_tic_n;
    irf_tk_tk_p = irf_tk_tk_n;
    irf_tl_tl_p = irf_tl_tl_n;
    irf_y_tc_p = irf_y_tc_n;
    irf_y_td_p = irf_y_td_n;
    irf_y_ti_p = irf_y_ti_n;
    irf_y_tic_p = irf_y_tic_n;
    irf_y_tk_p = irf_y_tk_n;
    irf_y_tl_p = irf_y_tl_n;
end


%% =========================================================================
%% PART 5: CALCULATE PRESENT VALUE MULTIPLIERS (Zubairy 2014 Style)
%% =========================================================================
% 
% Multiplier Definition (Zubairy 2014, p. 184):
%   PV Multiplier at horizon k = PV(DeltaY) / PV(DeltaF)
%
% where:
%   PV(DeltaY) = sum_{j=0}^{k-1} (1+R)^{-j} * DeltaY_{t+j}
%   PV(DeltaF) = sum_{j=0}^{k-1} (1+R)^{-j} * DeltaF_{t+j}
%
% For tax multipliers, DeltaF = change in tax REVENUE (not rate)
%
% NORMALIZATION: We express multipliers in terms of:
%   "$ change in GDP per $ change in fiscal variable"
%
% This matches Zubairy's "1% of GDP" normalization.
%% =========================================================================

fprintf('Calculating present-value multipliers...\n');

% Pre-allocate
mult_tau_c = zeros(n_horizons, 1);
mult_tau_l = zeros(n_horizons, 1);
mult_tau_k = zeros(n_horizons, 1);
mult_tau_i = zeros(n_horizons, 1);
mult_tau_d = zeros(n_horizons, 1);
mult_tau_ic = zeros(n_horizons, 1);

for h = 1:n_horizons
    k = horizons(h);
    discount_factors = (1+R).^(-(0:k-1)');
    
    % =============================================================
    % NOTE ON SIGN CONVENTIONS
    % For tax rates (tau_c, tau_l, tau_k, tau_i, tau_d):
    %   We report multipliers for TAX CUTS, so we flip the sign of all IRFs.
    %   The denominator is the present value of the REVENUE DECREASE, so we use -PV(DeltaRevenue).
    % For the investment tax credit (tau_ic):
    %   We report multipliers for an INCREASE in the credit (a larger subsidy).
    % =============================================================
    
    
    % =============================================================
    % Multipliers:
    %   For tax shocks: PV(DeltaY) / PV(Delta total tax revenue), tax cut convention
    %   For ITC: PV(DeltaY) / PV(Delta outlays), ITC increase convention
    % =============================================================

    % ===== CONSUMPTION TAX (CUT) =====
    y_path = y_ss + irf_y_tc_n(1:k)';
    rev_ss = Rev_total_ss;
    rev_path = rev_ss + irf_rev_total_tc_n(1:k)';
    dY = y_path - y_ss;
    dRev = rev_path - rev_ss;
    PV_Y = sum(discount_factors .* dY);
    PV_Rev = sum(discount_factors .* dRev);
    mult_tau_c(h) = PV_Y / (-PV_Rev);

    % ===== LABOR INCOME TAX (CUT) =====
    y_path = y_ss + irf_y_tl_n(1:k)';
    rev_ss = Rev_total_ss;
    rev_path = rev_ss + irf_rev_total_tl_n(1:k)';
    dY = y_path - y_ss;
    dRev = rev_path - rev_ss;
    PV_Y = sum(discount_factors .* dY);
    PV_Rev = sum(discount_factors .* dRev);
    mult_tau_l(h) = PV_Y / (-PV_Rev);

    % ===== CAPITAL INCOME TAX (CUT) =====
    y_path = y_ss + irf_y_tk_n(1:k)';
    rev_ss = Rev_total_ss;
    rev_path = rev_ss + irf_rev_total_tk_n(1:k)';
    dY = y_path - y_ss;
    dRev = rev_path - rev_ss;
    PV_Y = sum(discount_factors .* dY);
    PV_Rev = sum(discount_factors .* dRev);
    mult_tau_k(h) = PV_Y / (-PV_Rev);

    % ===== INTEREST INCOME TAX (CUT) =====
    y_path = y_ss + irf_y_ti_n(1:k)';
    rev_ss = Rev_total_ss;
    rev_path = rev_ss + irf_rev_total_ti_n(1:k)';
    dY = y_path - y_ss;
    dRev = rev_path - rev_ss;
    PV_Y = sum(discount_factors .* dY);
    PV_Rev = sum(discount_factors .* dRev);
    mult_tau_i(h) = PV_Y / (-PV_Rev);

    % ===== DIVIDEND TAX (CUT) =====
    y_path = y_ss + irf_y_td_n(1:k)';
    rev_ss = Rev_total_ss;
    rev_path = rev_ss + irf_rev_total_td_n(1:k)';
    dY = y_path - y_ss;
    dRev = rev_path - rev_ss;
    PV_Y = sum(discount_factors .* dY);
    PV_Rev = sum(discount_factors .* dRev);
    mult_tau_d(h) = PV_Y / (-PV_Rev);

    % ===== INVESTMENT TAX CREDIT (INCREASE) =====
    y_path = y_ss + irf_y_tic_n(1:k)';
    spend_ss = Spend_ic_ss;
    spend_path = spend_ss + irf_spend_ic_tic_n(1:k)';
    dY = y_path - y_ss;
    dSpend = spend_path - spend_ss;
    PV_Y = sum(discount_factors .* dY);
    PV_Spend = sum(discount_factors .* dSpend);
    mult_tau_ic(h) = PV_Y / PV_Spend;
end

fprintf('  Multipliers calculated for horizons: ');
fprintf('%d ', horizons);
fprintf('quarters\n\n');

%% =========================================================================
%% PART 6: DISPLAY MULTIPLIER TABLE (Zubairy 2014 Style)
%% =========================================================================

fprintf('=====================================================================\n');
fprintf('       PRESENT VALUE TAX MULTIPLIERS (Zubairy 2014 Style)\n');
fprintf('=====================================================================\n');
fprintf('\n');
fprintf('                            Quarter 1   Quarter 4   Quarter 12  Quarter 20\n');
fprintf('-------------------------------------------------------------------------\n');
fprintf('Consumption Tax (Cut)\n');
fprintf('  PV(DeltaY)/PV(DeltaT_total)  %8.2f    %8.2f    %8.2f    %8.2f\n', ...
    mult_tau_c(1), mult_tau_c(2), mult_tau_c(3), mult_tau_c(4));
fprintf('\n');
fprintf('Labor Tax (Cut)\n');
fprintf('  PV(DeltaY)/PV(DeltaT_total)  %8.2f    %8.2f    %8.2f    %8.2f\n', ...
    mult_tau_l(1), mult_tau_l(2), mult_tau_l(3), mult_tau_l(4));
fprintf('\n');
fprintf('Capital Tax (Cut)\n');
fprintf('  PV(DeltaY)/PV(DeltaT_total)  %8.2f    %8.2f    %8.2f    %8.2f\n', ...
    mult_tau_k(1), mult_tau_k(2), mult_tau_k(3), mult_tau_k(4));
fprintf('\n');
fprintf('Interest Tax (Cut)\n');
fprintf('  PV(DeltaY)/PV(DeltaT_total)  %8.2f    %8.2f    %8.2f    %8.2f\n', ...
    mult_tau_i(1), mult_tau_i(2), mult_tau_i(3), mult_tau_i(4));
fprintf('\n');
fprintf('Dividend Tax (Cut)\n');
fprintf('  PV(DeltaY)/PV(DeltaT_total)  %8.2f    %8.2f    %8.2f    %8.2f\n', ...
    mult_tau_d(1), mult_tau_d(2), mult_tau_d(3), mult_tau_d(4));
fprintf('\n');
fprintf('Investment Tax Credit (Increase)\n');
fprintf('  PV(DeltaY)/PV(DeltaOutlay) %8.2f    %8.2f    %8.2f    %8.2f\n', ...
    mult_tau_ic(1), mult_tau_ic(2), mult_tau_ic(3), mult_tau_ic(4));
fprintf('=========================================================================\n');
fprintf('\n');

% Benchmark comparison
fprintf('BENCHMARK COMPARISON:\n');
fprintf('-------------------------------------------------------------------------\n');
fprintf('Zubairy (2014) estimates:\n');
fprintf('  Labor Tax:               0.13        0.32        0.68        0.85\n');
fprintf('  Capital Tax:             0.34        0.43        0.52        0.46\n');
fprintf('\n');
fprintf('Mertens & Ravn (2013) estimates:\n');
fprintf('  Personal Income Tax:     Impact ~2.0, Peak ~2.5 at Q3\n');
fprintf('  Corporate Income Tax:    Impact ~0.4, Peak ~0.6 at Q4\n');
fprintf('=========================================================================\n');
fprintf('\n');

%% =========================================================================
%% PART 7: FIGURE 1 - TAX CUTS COMPARISON (Main Figure)
%% =========================================================================

figure('Position', [100 100 1400 900], 'Color', 'w', 'Name', 'Tax Cuts Comparison');

% Colors for different taxes
colors = struct();
colors.tc = [0.8 0.2 0.2];   % Red - Consumption tax
colors.tl = [0.2 0.2 0.8];   % Blue - Labor tax
colors.tk = [0.2 0.7 0.2];   % Green - Capital tax
colors.ti = [0.8 0.5 0.0];   % Orange - Interest tax
colors.td = [0.5 0.0 0.8];   % Purple - Dividend tax
colors.tic = [0.0 0.7 0.7];  % Cyan - ITC

% Output responses to tax CUTS (flip signs except ITC)
subplot(2,3,1);
plot(quarters, irf_y_tc_p*100, '-', 'Color', colors.tc, 'LineWidth', 2); hold on;
plot(quarters, irf_y_tl_p*100, '-', 'Color', colors.tl, 'LineWidth', 2);
plot(quarters, irf_y_tk_p*100, '-', 'Color', colors.tk, 'LineWidth', 2);
plot(quarters, irf_y_tic_p*100, '-', 'Color', colors.tic, 'LineWidth', 2);
plot(quarters, zeros(T,1), 'k--', 'LineWidth', 0.5);
xlabel('Quarters');
ylabel('Percent');
title('Output');
legend('\tau_c cut', '\tau_l cut', '\tau_k cut', 'ITC increase', 'Location', 'best', 'FontSize', 8);
xlim([0 T-1]);
grid on;

% Consumption responses
subplot(2,3,2);
plot(quarters, irf_c_tc_p*100, '-', 'Color', colors.tc, 'LineWidth', 2); hold on;
plot(quarters, irf_c_tl_p*100, '-', 'Color', colors.tl, 'LineWidth', 2);
plot(quarters, irf_c_tk_p*100, '-', 'Color', colors.tk, 'LineWidth', 2);
plot(quarters, irf_c_tic_p*100, '-', 'Color', colors.tic, 'LineWidth', 2);
plot(quarters, zeros(T,1), 'k--', 'LineWidth', 0.5);
xlabel('Quarters');
ylabel('Percent');
title('Consumption');
xlim([0 T-1]);
grid on;

% Investment responses
subplot(2,3,3);
plot(quarters, irf_i_tc_p*100, '-', 'Color', colors.tc, 'LineWidth', 2); hold on;
plot(quarters, irf_i_tl_p*100, '-', 'Color', colors.tl, 'LineWidth', 2);
plot(quarters, irf_i_tk_p*100, '-', 'Color', colors.tk, 'LineWidth', 2);
plot(quarters, irf_i_tic_p*100, '-', 'Color', colors.tic, 'LineWidth', 2);
plot(quarters, zeros(T,1), 'k--', 'LineWidth', 0.5);
xlabel('Quarters');
ylabel('Percent');
title('Investment');
xlim([0 T-1]);
grid on;

% Hours responses
subplot(2,3,4);
plot(quarters, irf_l_tc_p*100, '-', 'Color', colors.tc, 'LineWidth', 2); hold on;
plot(quarters, irf_l_tl_p*100, '-', 'Color', colors.tl, 'LineWidth', 2);
plot(quarters, irf_l_tk_p*100, '-', 'Color', colors.tk, 'LineWidth', 2);
plot(quarters, irf_l_tic_p*100, '-', 'Color', colors.tic, 'LineWidth', 2);
plot(quarters, zeros(T,1), 'k--', 'LineWidth', 0.5);
xlabel('Quarters');
ylabel('Percent');
title('Hours');
xlim([0 T-1]);
grid on;

% Interest and Dividend tax effects on output
subplot(2,3,5);
plot(quarters, irf_y_ti_p*100, '-', 'Color', colors.ti, 'LineWidth', 2); hold on;
plot(quarters, irf_y_td_p*100, '-', 'Color', colors.td, 'LineWidth', 2);
plot(quarters, zeros(T,1), 'k--', 'LineWidth', 0.5);
xlabel('Quarters');
ylabel('Percent');
title('Output: Interest & Dividend Tax Cuts');
legend('\tau_i cut', '\tau_d cut', 'Location', 'best', 'FontSize', 9);
xlim([0 T-1]);
grid on;

% Tax rate paths
subplot(2,3,6);
plot(quarters, irf_tc_tc_p*100, '-', 'Color', colors.tc, 'LineWidth', 2); hold on;
plot(quarters, irf_tl_tl_p*100, '-', 'Color', colors.tl, 'LineWidth', 2);
plot(quarters, irf_tk_tk_p*100, '-', 'Color', colors.tk, 'LineWidth', 2);
plot(quarters, zeros(T,1), 'k--', 'LineWidth', 0.5);
xlabel('Quarters');
ylabel('Percentage Points');
title('Tax Rate Paths (tax cuts, raw shock size)');
legend('\tau_c', '\tau_l', '\tau_k', 'Location', 'best', 'FontSize', 9);
xlim([0 T-1]);
grid on;

sgtitle('Impulse Responses to Tax Policy Shocks', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, 'IRF_tax_comparison.png');
saveas(gcf, 'IRF_tax_comparison.eps', 'epsc');

%% =========================================================================
%% PART 8: FIGURE 2 - CONSUMPTION TAX FOCUS (Key Result)
%% =========================================================================

figure('Position', [100 100 1000 700], 'Color', 'w', 'Name', 'Consumption Tax Shock');

subplot(2,2,1);
plot(quarters, irf_y_tc_p*100, 'r-', 'LineWidth', 2.5);
hold on;
plot(quarters, zeros(T,1), 'k--', 'LineWidth', 0.5);
xlabel('Quarters', 'FontSize', 11);
ylabel('Percent', 'FontSize', 11);
title('Output', 'FontSize', 12, 'FontWeight', 'bold');
xlim([0 T-1]);
grid on;
box on;

subplot(2,2,2);
plot(quarters, irf_c_tc_p*100, 'r-', 'LineWidth', 2.5);
hold on;
plot(quarters, zeros(T,1), 'k--', 'LineWidth', 0.5);
xlabel('Quarters', 'FontSize', 11);
ylabel('Percent', 'FontSize', 11);
title('Consumption', 'FontSize', 12, 'FontWeight', 'bold');
xlim([0 T-1]);
grid on;
box on;

subplot(2,2,3);
plot(quarters, irf_i_tc_p*100, 'r-', 'LineWidth', 2.5);
hold on;
plot(quarters, zeros(T,1), 'k--', 'LineWidth', 0.5);
xlabel('Quarters', 'FontSize', 11);
ylabel('Percent', 'FontSize', 11);
title('Investment', 'FontSize', 12, 'FontWeight', 'bold');
xlim([0 T-1]);
grid on;
box on;

subplot(2,2,4);
plot(quarters, irf_l_tc_p*100, 'r-', 'LineWidth', 2.5);
hold on;
plot(quarters, zeros(T,1), 'k--', 'LineWidth', 0.5);
xlabel('Quarters', 'FontSize', 11);
ylabel('Percent', 'FontSize', 11);
title('Hours', 'FontSize', 12, 'FontWeight', 'bold');
xlim([0 T-1]);
grid on;
box on;

sgtitle('Impulse Responses to Consumption Tax CUT', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, 'IRF_consumption_tax.png');
saveas(gcf, 'IRF_consumption_tax.eps', 'epsc');

%% =========================================================================
%% PART 9: FIGURE 3 - MULTIPLIERS OVER HORIZON
%% =========================================================================

figure('Position', [100 100 900 600], 'Color', 'w', 'Name', 'Multipliers Over Horizon');

% Calculate multipliers at each horizon
all_horizons = 1:T;
mult_tc_all = zeros(T, 1);
mult_tl_all = zeros(T, 1);
mult_tk_all = zeros(T, 1);
mult_ti_all = zeros(T, 1);
mult_td_all = zeros(T, 1);
mult_tic_all = zeros(T, 1);

for h = 1:T
    discount_factors = (1+R).^(-(0:h-1)');
    
    % Consumption tax (cut)
    sgn = -1;
    y_path = y_ss + sgn*irf_y_tc(1:h)';
    rev_ss = Rev_c_ss;
    rev_path = rev_ss + sgn*irf_rev_c_tc(1:h)';
    PV_Y = sum(discount_factors .* (y_path - y_ss));
    PV_Rev = sum(discount_factors .* (rev_path - rev_ss));
    mult_tc_all(h) = PV_Y / (-PV_Rev);
    
    % Labor income tax (cut)
    sgn = -1;
    y_path = y_ss + sgn*irf_y_tl(1:h)';
    p_path = p_ss + sgn*irf_p_tl(1:h)';
    w_path = w_ss + sgn*irf_w_tl(1:h)';
    l_path = l_ss + sgn*irf_l_tl(1:h)';
    tau_path = tau_l_ss + sgn*irf_tl_tl(1:h)';
    rev_ss = Rev_l_ss;
    rev_path = rev_ss + sgn*irf_rev_l_tl(1:h)';
    PV_Y = sum(discount_factors .* (y_path - y_ss));
    PV_Rev = sum(discount_factors .* (rev_path - rev_ss));
    mult_tl_all(h) = PV_Y / (-PV_Rev);
    
    % Capital income tax (cut)
    sgn = -1;
    y_path = y_ss + sgn*irf_y_tk(1:h)';
    rk_path = r_k_ss + sgn*irf_rk_tk(1:h)';
    v_path  = v_ss  + sgn*irf_v_tk(1:h)';
    k_path  = k_ss  + sgn*irf_k_tk(1:h)';
    ktau_path = k_tau_ss + sgn*irf_ktau_tk(1:h)';
    i_path  = i_ss  + sgn*irf_i_tk(1:h)';
    tau_path = tau_k_ss + sgn*irf_tk_tk(1:h)';
    rev_ss = Rev_k_ss;
    rev_path = rev_ss + sgn*irf_rev_k_tk(1:h)';
    PV_Y = sum(discount_factors .* (y_path - y_ss));
    PV_Rev = sum(discount_factors .* (rev_path - rev_ss));
    mult_tk_all(h) = PV_Y / (-PV_Rev);
    
    % Interest income tax (cut)
    sgn = -1;
    y_path = y_ss + sgn*irf_y_ti(1:h)';
    rev_ss = Rev_i_ss;
    rev_path = rev_ss + sgn*irf_rev_i_ti(1:h)';
    PV_Y = sum(discount_factors .* (y_path - y_ss));
    PV_Rev = sum(discount_factors .* (rev_path - rev_ss));
    mult_ti_all(h) = PV_Y / (-PV_Rev);

    % Dividend tax (cut)
    sgn = -1;
    y_path = y_ss + sgn*irf_y_td(1:h)';
    rev_ss = Rev_d_ss;
    rev_path = rev_ss + sgn*irf_rev_d_td(1:h)';
    PV_Y = sum(discount_factors .* (y_path - y_ss));
    PV_Rev = sum(discount_factors .* (rev_path - rev_ss));
    mult_td_all(h) = PV_Y / (-PV_Rev);

    % ITC (increase)
    sgn = 1;
    y_path = y_ss + sgn*irf_y_tic(1:h)';
    i_path = i_ss + sgn*irf_i_tic(1:h)';
    tau_path = tau_ic_ss + sgn*irf_tic_tic(1:h)';
    spend_ss = Spend_ic_ss;
    spend_path = spend_ss + sgn*irf_spend_ic_tic(1:h)';
    PV_Y = sum(discount_factors .* (y_path - y_ss));
    PV_Spend = sum(discount_factors .* (spend_path - spend_ss));
    mult_tic_all(h) = PV_Y / PV_Spend;
end

plot(all_horizons, mult_tc_all, '-', 'Color', colors.tc, 'LineWidth', 2.5); hold on;
plot(all_horizons, mult_tl_all, '-', 'Color', colors.tl, 'LineWidth', 2);
plot(all_horizons, mult_tk_all, '-', 'Color', colors.tk, 'LineWidth', 2);
plot(all_horizons, mult_ti_all, '-', 'Color', colors.ti, 'LineWidth', 2);
plot(all_horizons, mult_td_all, '-', 'Color', colors.td, 'LineWidth', 2);
plot(all_horizons, mult_tic_all, '-', 'Color', colors.tic, 'LineWidth', 2);
plot(all_horizons, ones(T,1), 'k--', 'LineWidth', 0.5);

xlabel('Horizon (Quarters)', 'FontSize', 12);
ylabel('Present Value Multiplier', 'FontSize', 12);
title('Tax Multipliers by Horizon', 'FontSize', 14, 'FontWeight', 'bold');
legend('Cons. Tax Cut', 'Labor Tax Cut', 'Capital Tax Cut', 'Interest Tax Cut', 'Dividend Tax Cut', 'ITC Increase', ...
    'Location', 'best', 'FontSize', 10);
xlim([1 T]);
grid on;
box on;

saveas(gcf, 'multipliers_by_horizon.png');
saveas(gcf, 'multipliers_by_horizon.eps', 'epsc');

%% =========================================================================
%% PART 10: SAVE RESULTS
%% =========================================================================

% Save multiplier table to CSV
T_mult = table(horizons', mult_tau_c, mult_tau_l, mult_tau_k, ...
    mult_tau_i, mult_tau_d, mult_tau_ic, ...
    'VariableNames', {'Horizon_Q', 'ConsTax', 'LaborTax', ...
    'CapitalTax', 'InterestTax', 'DividendTax', 'InvTaxCredit'});
writetable(T_mult, 'multiplier_table.csv');

% Save all results to MAT file
save('fiscal_multiplier_results.mat', ...
    'irf_*', 'mult_*', 'horizons', 'quarters', '*_ss', 'T_*_base');

fprintf('\n');
fprintf('Results saved to:\n');
fprintf('  - multiplier_table.csv\n');
fprintf('  - fiscal_multiplier_results.mat\n');
fprintf('  - IRF_tax_comparison.png/eps\n');
fprintf('  - IRF_consumption_tax.png/eps\n');
fprintf('  - multipliers_by_horizon.png/eps\n');
fprintf('\n');
fprintf('Done!\n');
fprintf('================================================================\n');