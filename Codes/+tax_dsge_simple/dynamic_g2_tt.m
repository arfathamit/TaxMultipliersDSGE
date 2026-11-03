function T = dynamic_g2_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_g2_tt(T, y, x, params, steady_state, it_)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double  vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double  vector of endogenous variables in the order stored
%                                                    in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double  matrix of exogenous variables (in declaration order)
%                                                    for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double  vector of steady state values
%   params        [M_.param_nbr by 1]        double  vector of parameter values in declaration order
%   it_           scalar                     double  time period for exogenous variables for which
%                                                    to evaluate the model
%
% Output:
%   T           [#temp variables by 1]       double  vector of temporary terms
%

assert(length(T) >= 79);

T = tax_dsge_simple.dynamic_g1_tt(T, y, x, params, steady_state, it_);

T(72) = (-((-y(16))*(y(4)+y(4))))/(y(4)*y(4)*y(4)*y(4));
T(73) = getPowerDeriv(y(13)/T(17),1/(params(11)-1),2);
T(74) = T(47)*T(73);
T(75) = getPowerDeriv(y(14),1-params(11),2);
T(76) = getPowerDeriv(y(9)*y(3),params(11),2);
T(77) = getPowerDeriv(T(30),(-y(26)),2);
T(78) = exp(y(25))^(T(27)-1);
T(79) = T(30)^((-y(26))-1);

end
