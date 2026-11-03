%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

if isoctave || matlab_ver_less_than('8.6')
    clear all
else
    clearvars -global
    clear_persistent_variables(fileparts(which('dynare')), false)
end
tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'tax_dsge';
M_.dynare_version = '5.4';
oo_.dynare_version = '5.4';
options_.dynare_version = '5.4';
%
% Some global variables initialization
%
global_initialization;
M_.exo_names = cell(3,1);
M_.exo_names_tex = cell(3,1);
M_.exo_names_long = cell(3,1);
M_.exo_names(1) = {'e_b'};
M_.exo_names_tex(1) = {'e\_b'};
M_.exo_names_long(1) = {'e_b'};
M_.exo_names(2) = {'e_p'};
M_.exo_names_tex(2) = {'e\_p'};
M_.exo_names_long(2) = {'e_p'};
M_.exo_names(3) = {'e_z'};
M_.exo_names_tex(3) = {'e\_z'};
M_.exo_names_long(3) = {'e_z'};
M_.endo_names = cell(22,1);
M_.endo_names_tex = cell(22,1);
M_.endo_names_long = cell(22,1);
M_.endo_names(1) = {'p'};
M_.endo_names_tex(1) = {'p'};
M_.endo_names_long(1) = {'p'};
M_.endo_names(2) = {'p_star'};
M_.endo_names_tex(2) = {'p\_star'};
M_.endo_names_long(2) = {'p_star'};
M_.endo_names(3) = {'v'};
M_.endo_names_tex(3) = {'v'};
M_.endo_names_long(3) = {'v'};
M_.endo_names(4) = {'z'};
M_.endo_names_tex(4) = {'z'};
M_.endo_names_long(4) = {'z'};
M_.endo_names(5) = {'epsilon_b'};
M_.endo_names_tex(5) = {'epsilon\_b'};
M_.endo_names_long(5) = {'epsilon_b'};
M_.endo_names(6) = {'lambda_p'};
M_.endo_names_tex(6) = {'lambda\_p'};
M_.endo_names_long(6) = {'lambda_p'};
M_.endo_names(7) = {'r'};
M_.endo_names_tex(7) = {'r'};
M_.endo_names_long(7) = {'r'};
M_.endo_names(8) = {'r_k'};
M_.endo_names_tex(8) = {'r\_k'};
M_.endo_names_long(8) = {'r_k'};
M_.endo_names(9) = {'w'};
M_.endo_names_tex(9) = {'w'};
M_.endo_names_long(9) = {'w'};
M_.endo_names(10) = {'l'};
M_.endo_names_tex(10) = {'l'};
M_.endo_names_long(10) = {'l'};
M_.endo_names(11) = {'k'};
M_.endo_names_tex(11) = {'k'};
M_.endo_names_long(11) = {'k'};
M_.endo_names(12) = {'k_tau'};
M_.endo_names_tex(12) = {'k\_tau'};
M_.endo_names_long(12) = {'k_tau'};
M_.endo_names(13) = {'i'};
M_.endo_names_tex(13) = {'i'};
M_.endo_names_long(13) = {'i'};
M_.endo_names(14) = {'y'};
M_.endo_names_tex(14) = {'y'};
M_.endo_names_long(14) = {'y'};
M_.endo_names(15) = {'d'};
M_.endo_names_tex(15) = {'d'};
M_.endo_names_long(15) = {'d'};
M_.endo_names(16) = {'c'};
M_.endo_names_tex(16) = {'c'};
M_.endo_names_long(16) = {'c'};
M_.endo_names(17) = {'lambda'};
M_.endo_names_tex(17) = {'lambda'};
M_.endo_names_long(17) = {'lambda'};
M_.endo_names(18) = {'q'};
M_.endo_names_tex(18) = {'q'};
M_.endo_names_long(18) = {'q'};
M_.endo_names(19) = {'g'};
M_.endo_names_tex(19) = {'g'};
M_.endo_names_long(19) = {'g'};
M_.endo_names(20) = {'b'};
M_.endo_names_tex(20) = {'b'};
M_.endo_names_long(20) = {'b'};
M_.endo_names(21) = {'x'};
M_.endo_names_tex(21) = {'x'};
M_.endo_names_long(21) = {'x'};
M_.endo_names(22) = {'mc'};
M_.endo_names_tex(22) = {'mc'};
M_.endo_names_long(22) = {'mc'};
M_.endo_partitions = struct();
M_.param_names = cell(36,1);
M_.param_names_tex = cell(36,1);
M_.param_names_long = cell(36,1);
M_.param_names(1) = {'betta'};
M_.param_names_tex(1) = {'betta'};
M_.param_names_long(1) = {'betta'};
M_.param_names(2) = {'siggma'};
M_.param_names_tex(2) = {'siggma'};
M_.param_names_long(2) = {'siggma'};
M_.param_names(3) = {'zetta'};
M_.param_names_tex(3) = {'zetta'};
M_.param_names_long(3) = {'zetta'};
M_.param_names(4) = {'siggma_g'};
M_.param_names_tex(4) = {'siggma\_g'};
M_.param_names_long(4) = {'siggma_g'};
M_.param_names(5) = {'chi_g'};
M_.param_names_tex(5) = {'chi\_g'};
M_.param_names_long(5) = {'chi_g'};
M_.param_names(6) = {'rho_b'};
M_.param_names_tex(6) = {'rho\_b'};
M_.param_names_long(6) = {'rho_b'};
M_.param_names(7) = {'gamma'};
M_.param_names_tex(7) = {'gamma'};
M_.param_names_long(7) = {'gamma'};
M_.param_names(8) = {'delta_0'};
M_.param_names_tex(8) = {'delta\_0'};
M_.param_names_long(8) = {'delta_0'};
M_.param_names(9) = {'delta_1'};
M_.param_names_tex(9) = {'delta\_1'};
M_.param_names_long(9) = {'delta_1'};
M_.param_names(10) = {'delta_2'};
M_.param_names_tex(10) = {'delta\_2'};
M_.param_names_long(10) = {'delta_2'};
M_.param_names(11) = {'alfa'};
M_.param_names_tex(11) = {'alfa'};
M_.param_names_long(11) = {'alfa'};
M_.param_names(12) = {'lambda_bar_p'};
M_.param_names_tex(12) = {'lambda\_bar\_p'};
M_.param_names_long(12) = {'lambda_bar_p'};
M_.param_names(13) = {'theta'};
M_.param_names_tex(13) = {'theta'};
M_.param_names_long(13) = {'theta'};
M_.param_names(14) = {'phi_1'};
M_.param_names_tex(14) = {'phi\_1'};
M_.param_names_long(14) = {'phi_1'};
M_.param_names(15) = {'phi_2'};
M_.param_names_tex(15) = {'phi\_2'};
M_.param_names_long(15) = {'phi_2'};
M_.param_names(16) = {'rho_r'};
M_.param_names_tex(16) = {'rho\_r'};
M_.param_names_long(16) = {'rho_r'};
M_.param_names(17) = {'rho_z'};
M_.param_names_tex(17) = {'rho\_z'};
M_.param_names_long(17) = {'rho_z'};
M_.param_names(18) = {'rho_lambda'};
M_.param_names_tex(18) = {'rho\_lambda'};
M_.param_names_long(18) = {'rho_lambda'};
M_.param_names(19) = {'siggma_z'};
M_.param_names_tex(19) = {'siggma\_z'};
M_.param_names_long(19) = {'siggma_z'};
M_.param_names(20) = {'siggma_b'};
M_.param_names_tex(20) = {'siggma\_b'};
M_.param_names_long(20) = {'siggma_b'};
M_.param_names(21) = {'siggma_p'};
M_.param_names_tex(21) = {'siggma\_p'};
M_.param_names_long(21) = {'siggma_p'};
M_.param_names(22) = {'tau_c'};
M_.param_names_tex(22) = {'tau\_c'};
M_.param_names_long(22) = {'tau_c'};
M_.param_names(23) = {'tau_i'};
M_.param_names_tex(23) = {'tau\_i'};
M_.param_names_long(23) = {'tau_i'};
M_.param_names(24) = {'tau_k'};
M_.param_names_tex(24) = {'tau\_k'};
M_.param_names_long(24) = {'tau_k'};
M_.param_names(25) = {'tau_l'};
M_.param_names_tex(25) = {'tau\_l'};
M_.param_names_long(25) = {'tau_l'};
M_.param_names(26) = {'tau_d'};
M_.param_names_tex(26) = {'tau\_d'};
M_.param_names_long(26) = {'tau_d'};
M_.param_names(27) = {'tau_ic'};
M_.param_names_tex(27) = {'tau\_ic'};
M_.param_names_long(27) = {'tau_ic'};
M_.param_names(28) = {'e_tau'};
M_.param_names_tex(28) = {'e\_tau'};
M_.param_names_long(28) = {'e_tau'};
M_.param_names(29) = {'delta_tau'};
M_.param_names_tex(29) = {'delta\_tau'};
M_.param_names_long(29) = {'delta_tau'};
M_.param_names(30) = {'gamma_x'};
M_.param_names_tex(30) = {'gamma\_x'};
M_.param_names_long(30) = {'gamma_x'};
M_.param_names(31) = {'lbar'};
M_.param_names_tex(31) = {'lbar'};
M_.param_names_long(31) = {'lbar'};
M_.param_names(32) = {'ybar'};
M_.param_names_tex(32) = {'ybar'};
M_.param_names_long(32) = {'ybar'};
M_.param_names(33) = {'rbar'};
M_.param_names_tex(33) = {'rbar'};
M_.param_names_long(33) = {'rbar'};
M_.param_names(34) = {'r_k_bar'};
M_.param_names_tex(34) = {'r\_k\_bar'};
M_.param_names_long(34) = {'r_k_bar'};
M_.param_names(35) = {'wbar'};
M_.param_names_tex(35) = {'wbar'};
M_.param_names_long(35) = {'wbar'};
M_.param_names(36) = {'mc_bar'};
M_.param_names_tex(36) = {'mc\_bar'};
M_.param_names_long(36) = {'mc_bar'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 3;
M_.endo_nbr = 22;
M_.param_nbr = 36;
M_.orig_endo_nbr = 22;
M_.aux_vars = [];
M_ = setup_solvers(M_);
M_.Sigma_e = zeros(3, 3);
M_.Correlation_matrix = eye(3, 3);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = true;
M_.det_shocks = [];
M_.surprise_shocks = [];
M_.heteroskedastic_shocks.Qvalue_orig = [];
M_.heteroskedastic_shocks.Qscale_orig = [];
options_.linear = false;
options_.block = false;
options_.bytecode = false;
options_.use_dll = false;
M_.nonzero_hessian_eqs = [1 2 3 4 5 6 7 8 11 12 13 14 15 16 18 19 21];
M_.hessian_eq_zero = isempty(M_.nonzero_hessian_eqs);
M_.orig_eq_nbr = 22;
M_.eq_nbr = 22;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./+' M_.fname '/set_auxiliary_variables.m'], 'file') == 2;
M_.epilogue_names = {};
M_.epilogue_var_list_ = {};
M_.orig_maximum_endo_lag = 1;
M_.orig_maximum_endo_lead = 1;
M_.orig_maximum_exo_lag = 0;
M_.orig_maximum_exo_lead = 0;
M_.orig_maximum_exo_det_lag = 0;
M_.orig_maximum_exo_det_lead = 0;
M_.orig_maximum_lag = 1;
M_.orig_maximum_lead = 1;
M_.orig_maximum_lag_with_diffs_expanded = 1;
M_.lead_lag_incidence = [
 1 9 31;
 0 10 32;
 0 11 33;
 2 12 0;
 3 13 0;
 4 14 0;
 0 15 34;
 0 16 35;
 0 17 0;
 0 18 0;
 5 19 0;
 6 20 0;
 7 21 36;
 0 22 0;
 0 23 0;
 0 24 0;
 0 25 37;
 0 26 38;
 0 27 0;
 8 28 0;
 0 29 0;
 0 30 0;]';
M_.nstatic = 8;
M_.nfwrd   = 6;
M_.npred   = 6;
M_.nboth   = 2;
M_.nsfwrd   = 8;
M_.nspred   = 8;
M_.ndynamic   = 14;
M_.dynamic_tmp_nbr = [38; 38; 15; 0; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , 'c' ;
  2 , 'name' , 'lambda' ;
  3 , 'name' , '3' ;
  4 , 'name' , 'q' ;
  5 , 'name' , '5' ;
  6 , 'name' , '6' ;
  7 , 'name' , '7' ;
  8 , 'name' , 'k' ;
  9 , 'name' , 'k_tau' ;
  10 , 'name' , 'y' ;
  11 , 'name' , 'g' ;
  12 , 'name' , 'x' ;
  13 , 'name' , 'r' ;
  14 , 'name' , '14' ;
  15 , 'name' , 'l' ;
  16 , 'name' , 'mc' ;
  17 , 'name' , 'p_star' ;
  18 , 'name' , 'd' ;
  19 , 'name' , 'p' ;
  20 , 'name' , 'z' ;
  21 , 'name' , '21' ;
  22 , 'name' , 'epsilon_b' ;
};
M_.mapping.p.eqidx = [1 3 4 6 7 11 12 13 14 15 16 17 18 19 ];
M_.mapping.p_star.eqidx = [17 18 19 ];
M_.mapping.v.eqidx = [1 4 7 8 11 12 14 ];
M_.mapping.z.eqidx = [16 20 ];
M_.mapping.epsilon_b.eqidx = [2 3 22 ];
M_.mapping.lambda_p.eqidx = [17 18 19 21 ];
M_.mapping.r.eqidx = [1 6 11 12 13 ];
M_.mapping.r_k.eqidx = [1 4 7 11 12 14 16 ];
M_.mapping.w.eqidx = [1 3 11 12 15 16 ];
M_.mapping.l.eqidx = [1 2 3 11 12 15 ];
M_.mapping.k.eqidx = [1 7 8 11 12 14 ];
M_.mapping.k_tau.eqidx = [1 9 11 12 ];
M_.mapping.i.eqidx = [1 5 8 9 10 11 12 ];
M_.mapping.y.eqidx = [10 13 14 15 18 ];
M_.mapping.d.eqidx = [1 11 12 18 ];
M_.mapping.c.eqidx = [1 2 3 10 11 12 ];
M_.mapping.lambda.eqidx = [2 3 4 5 6 7 ];
M_.mapping.q.eqidx = [4 5 7 ];
M_.mapping.g.eqidx = [10 11 ];
M_.mapping.b.eqidx = [1 11 12 ];
M_.mapping.x.eqidx = [1 12 ];
M_.mapping.mc.eqidx = [16 17 18 ];
M_.mapping.e_b.eqidx = [22 ];
M_.mapping.e_p.eqidx = [21 ];
M_.mapping.e_z.eqidx = [20 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [1 4 5 6 11 12 13 20 ];
M_.exo_names_orig_ord = [1:3];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(22, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(3, 1);
M_.params = NaN(36, 1);
M_.endo_trends = struct('deflator', cell(22, 1), 'log_deflator', cell(22, 1), 'growth_factor', cell(22, 1), 'log_growth_factor', cell(22, 1));
M_.NNZDerivatives = [132; 258; -1; ];
M_.static_tmp_nbr = [22; 5; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
M_.params(1) = 0.95;
betta = M_.params(1);
M_.params(2) = 2;
siggma = M_.params(2);
M_.params(3) = 0.29;
zetta = M_.params(3);
M_.params(4) = 1.1;
siggma_g = M_.params(4);
M_.params(5) = 0.5;
chi_g = M_.params(5);
M_.params(6) = 0.9;
rho_b = M_.params(6);
M_.params(7) = 2;
gamma = M_.params(7);
M_.params(8) = 0.1;
delta_0 = M_.params(8);
M_.params(9) = 0.152632;
delta_1 = M_.params(9);
M_.params(10) = 0.01;
delta_2 = M_.params(10);
M_.params(11) = 0.33;
alfa = M_.params(11);
M_.params(12) = 8;
lambda_bar_p = M_.params(12);
M_.params(13) = 0.82;
theta = M_.params(13);
M_.params(17) = 0.7;
rho_z = M_.params(17);
M_.params(18) = 0;
rho_lambda = M_.params(18);
M_.params(19) = 0.01;
siggma_z = M_.params(19);
M_.params(21) = 0.01;
siggma_p = M_.params(21);
M_.params(20) = 0.01;
siggma_b = M_.params(20);
M_.params(14) = 0.2;
phi_1 = M_.params(14);
M_.params(15) = 0.03;
phi_2 = M_.params(15);
M_.params(16) = 0.1;
rho_r = M_.params(16);
M_.params(22) = 0.1;
tau_c = M_.params(22);
M_.params(23) = 0.2;
tau_i = M_.params(23);
M_.params(24) = 0.35;
tau_k = M_.params(24);
M_.params(25) = 0.25;
tau_l = M_.params(25);
M_.params(26) = 0.15;
tau_d = M_.params(26);
M_.params(27) = 0;
tau_ic = M_.params(27);
M_.params(28) = 0;
e_tau = M_.params(28);
M_.params(29) = 0.1;
delta_tau = M_.params(29);
M_.params(30) = 0.2;
gamma_x = M_.params(30);
M_.params(31) = 0.5;
lbar = M_.params(31);
M_.params(33) = (1-M_.params(1))/(M_.params(1)*(1-M_.params(23)));
rbar = M_.params(33);
M_.params(34) = (1/M_.params(1)-1+M_.params(8))/(1-M_.params(24))*(1-M_.params(27)-M_.params(24)*M_.params(28)-M_.params(1)*M_.params(24)*M_.params(29)*(1-M_.params(28)));
r_k_bar = M_.params(34);
M_.params(35) = (1-M_.params(11))*(M_.params(34)/M_.params(11))^((-M_.params(11))/(1-M_.params(11)));
wbar = M_.params(35);
M_.params(32) = (M_.params(11)/(1-M_.params(11))*M_.params(35)/M_.params(34))^M_.params(11)*M_.params(31);
ybar = M_.params(32);
M_.params(36) = M_.params(35)^(1-M_.params(11))*M_.params(34)^M_.params(11)/(M_.params(11)^M_.params(11)*(1-M_.params(11))^(1-M_.params(11)));
mc_bar = M_.params(36);
%
% INITVAL instructions
%
options_.initval_file = false;
oo_.steady_state(1) = 1;
oo_.steady_state(2) = 1;
oo_.steady_state(3) = 1;
oo_.steady_state(4) = 0;
oo_.steady_state(5) = 0;
oo_.steady_state(6) = 1;
oo_.steady_state(7) = 0.06;
oo_.steady_state(8) = 0.22;
oo_.steady_state(9) = 0.1;
oo_.steady_state(10) = 0.3;
oo_.steady_state(11) = 1;
oo_.steady_state(12) = 1;
oo_.steady_state(13) = 0.1;
oo_.steady_state(14) = 0.446;
oo_.steady_state(15) = 0;
oo_.steady_state(16) = .246;
oo_.steady_state(17) = 1;
oo_.steady_state(18) = 0.95;
oo_.steady_state(19) = 0.2;
oo_.steady_state(20) = 0.1;
oo_.steady_state(21) = 0.1;
oo_.steady_state(22) = 1;
oo_.exo_steady_state(1) = 0;
oo_.exo_steady_state(2) = 0;
oo_.exo_steady_state(3) = 0;
if M_.exo_nbr > 0
	oo_.exo_simul = ones(M_.maximum_lag,1)*oo_.exo_steady_state';
end
if M_.exo_det_nbr > 0
	oo_.exo_det_simul = ones(M_.maximum_lag,1)*oo_.exo_det_steady_state';
end
steady;
oo_.dr.eigval = check(M_,options_,oo_);
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = M_.params(20)^2;
M_.Sigma_e(2, 2) = M_.params(21)^2;
M_.Sigma_e(3, 3) = M_.params(19)^2;
options_.order = 2;
var_list_ = {};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);


oo_.time = toc(tic0);
disp(['Total computing time : ' dynsec2hms(oo_.time) ]);
if ~exist([M_.dname filesep 'Output'],'dir')
    mkdir(M_.dname,'Output');
end
save([M_.dname filesep 'Output' filesep 'tax_dsge_results.mat'], 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'tax_dsge_results.mat'], 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'tax_dsge_results.mat'], 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'tax_dsge_results.mat'], 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'tax_dsge_results.mat'], 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'tax_dsge_results.mat'], 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'tax_dsge_results.mat'], 'oo_recursive_', '-append');
end
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
