function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
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

assert(length(T) >= 71);

T = tax_dsge_simple.dynamic_resid_tt(T, y, x, params, steady_state, it_);

T(32) = (-(y(13)*T(17)))/(T(17)*T(17));
T(33) = getPowerDeriv(y(13)/T(17),1/(params(11)-1),1);
T(34) = params(11)*exp(y(8))*getPowerDeriv(y(9),params(11),1);
T(35) = T(16)*T(34);
T(36) = (-(y(13)*T(35)));
T(37) = T(36)/(T(17)*T(17));
T(38) = getPowerDeriv(y(9)*y(3),params(11),1);
T(39) = 1/T(17);
T(40) = getPowerDeriv(y(13),params(11),1);
T(41) = T(20)*T(40);
T(42) = (-(getPowerDeriv(1-y(14),(1-params(2))*(1-params(3)),1)));
T(43) = (-(exp(y(10))*T(1)*T(42)));
T(44) = (-(getPowerDeriv(1-y(14),(1-params(2))*(1-params(3))-1,1)));
T(45) = exp(y(10))*T(4)*T(44);
T(46) = getPowerDeriv(y(14),1-params(11),1);
T(47) = (-(y(13)*T(15)*T(46)))/(T(17)*T(17));
T(48) = (-y(16))/(y(4)*y(4));
T(49) = (-(T(9)*T(48)*2*(T(10)-1)));
T(50) = 1/y(4);
T(51) = (-(T(9)*2*(T(10)-1)*T(50)));
T(52) = (-y(31))/(y(16)*y(16));
T(53) = 2*y(31)/y(16);
T(54) = T(52)*T(53);
T(55) = 1/y(16);
T(56) = T(53)*T(55);
T(57) = getPowerDeriv(y(18),1-params(11),1);
T(58) = T(21)*T(57);
T(59) = params(3)*getPowerDeriv(y(20),params(3)*(1-params(2))-1,1);
T(60) = (-(exp(y(10))*T(2)*T(59)));
T(61) = (1-params(3))*getPowerDeriv(y(20),params(3)*(1-params(2)),1);
T(62) = exp(y(10))*T(5)*T(61);
T(63) = exp(y(25))*getPowerDeriv(exp(y(25)),T(27),1);
T(64) = getPowerDeriv(T(30),(-y(26)),1);
T(65) = 1/(y(26)*y(26));
T(66) = T(65)*log(y(7));
T(67) = T(65)*log(exp(y(25)));
T(68) = params(13)*T(28)*T(66)+(1-params(13))*T(29)*T(67);
T(69) = (-y(26))*T(68)/T(30)-log(T(30));
T(70) = getPowerDeriv(y(7),T(27),1);
T(71) = params(13)*T(70);

end
