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

assert(length(T) >= 91);

T = tax_dsge.dynamic_g1_tt(T, y, x, params, steady_state, it_);

T(77) = (-((-((1+y(15)*(1-params(23)))*y(8)))*(y(9)+y(9))))/(y(9)*y(9)*y(9)*y(9));
T(78) = (-1)/(y(9)*y(9));
T(79) = (-((-y(21))*(y(7)+y(7))))/(y(7)*y(7)*y(7)*y(7));
T(80) = (-((-y(28))*(y(9)+y(9))))/(y(9)*y(9)*y(9)*y(9))+(-((-(y(17)*params(25)*y(18)))*(y(9)+y(9))))/(y(9)*y(9)*y(9)*y(9))+(-((-(params(26)*y(23)))*(y(9)+y(9))))/(y(9)*y(9)*y(9)*y(9))-T(77);
T(81) = getPowerDeriv(1+(y(9)-y(1))/y(1),params(14),2);
T(82) = getPowerDeriv(T(27),T(28),2);
T(83) = T(27)^(T(28)-1);
T(84) = getPowerDeriv(T(31),T(28),2);
T(85) = T(31)^(T(28)-1);
T(86) = (-(y(14)+y(14)))/(y(14)*y(14)*y(14)*y(14));
T(87) = getPowerDeriv(T(37),(-y(14)),2);
T(88) = y(1)^(T(34)-1);
T(89) = T(37)^((-y(14))-1);
T(90) = (-y(14))*(T(52)*((-y(14))-1)/T(37)-log(T(37)))*T(89)-T(89);
T(91) = y(10)^(T(34)-1);

end
