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

assert(length(T) >= 76);

T = tax_dsge.dynamic_resid_tt(T, y, x, params, steady_state, it_);

T(39) = ((-y(1))-(y(9)-y(1)))/(y(1)*y(1));
T(40) = getPowerDeriv(1+(y(9)-y(1))/y(1),params(14),1);
T(41) = 1/y(9);
T(42) = getPowerDeriv(T(27),T(28),1);
T(43) = params(13)*getPowerDeriv(y(1),T(34),1);
T(44) = getPowerDeriv(T(37),(-y(14)),1);
T(45) = (-y(28))/(y(9)*y(9))+(-(y(17)*params(25)*y(18)))/(y(9)*y(9))+(-(params(26)*y(23)))/(y(9)*y(9))-(-((1+y(15)*(1-params(23)))*y(8)))/(y(9)*y(9));
T(46) = T(26)*T(26);
T(47) = (-y(1))/(y(9)*y(9));
T(48) = (-y(10))/(y(9)*y(9));
T(49) = getPowerDeriv(T(31),T(28),1);
T(50) = (1-params(13))*getPowerDeriv(y(10),T(34),1);
T(51) = 1/(y(14)*y(14));
T(52) = params(13)*T(35)*T(51)*log(y(1))+(1-params(13))*T(36)*T(51)*log(y(10));
T(53) = (-y(14))*T(52)/T(37)-log(T(37));
T(54) = 1/(1+params(33))*getPowerDeriv((1+y(15))/(1+params(33)),params(16),1);
T(55) = getPowerDeriv(y(16),params(11),1);
T(56) = T(21)*T(55);
T(57) = getPowerDeriv(y(17),1-params(11),1);
T(58) = T(22)*T(57);
T(59) = (-(getPowerDeriv(1-y(18),(1-params(2))*(1-params(3)),1)));
T(60) = (-(exp(y(13))*T(1)*T(59)/(1+params(22))));
T(61) = (-(getPowerDeriv(1-y(18),(1-params(2))*(1-params(3))-1,1)));
T(62) = exp(y(13))*T(4)*T(61);
T(63) = (-y(21))/(y(7)*y(7));
T(64) = (-(T(9)*T(63)*2*(T(10)-1)));
T(65) = 1/y(7);
T(66) = (-(T(9)*2*(T(10)-1)*T(65)));
T(67) = (-y(36))/(y(21)*y(21));
T(68) = 2*y(36)/y(21);
T(69) = T(67)*T(68);
T(70) = 1/y(21);
T(71) = T(68)*T(70);
T(72) = 1/params(32)*getPowerDeriv(y(22)/params(32),params(15),1);
T(73) = params(3)*getPowerDeriv(y(24),params(3)*(1-params(2))-1,1);
T(74) = (-(exp(y(13))*T(2)*T(73)/(1+params(22))));
T(75) = (1-params(3))*getPowerDeriv(y(24),params(3)*(1-params(2)),1);
T(76) = exp(y(13))*T(5)*T(75);

end
