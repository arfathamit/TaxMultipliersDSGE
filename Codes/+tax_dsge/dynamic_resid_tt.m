function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
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

assert(length(T) >= 38);

T(1) = params(3)*y(24)^(params(3)*(1-params(2))-1);
T(2) = (1-y(18))^((1-params(2))*(1-params(3)));
T(3) = exp(y(13))*T(1)*T(2)/(1+params(22));
T(4) = (1-params(3))*y(24)^(params(3)*(1-params(2)));
T(5) = (1-y(18))^((1-params(2))*(1-params(3))-1);
T(6) = exp(y(13))*T(4)*T(5);
T(7) = params(10)/2;
T(8) = 1-params(8)-params(9)*(y(33)-1)-T(7)*(y(33)-1)^2;
T(9) = params(7)/2;
T(10) = y(21)/y(7);
T(11) = 1-T(9)*(T(10)-1)^2;
T(12) = y(38)*params(7)*(y(36)/y(21)-1);
T(13) = (y(36)/y(21))^2;
T(14) = 1-(params(8)+params(9)*(y(11)-1)+T(7)*(y(11)-1)^2);
T(15) = y(28)/y(9)+y(24)*params(22)+y(17)*params(25)*y(18)/y(9)+y(5)*y(11)*params(24)*y(16)+params(26)*y(23)/y(9)-params(27)*y(21)-params(24)*params(29)*y(6)-y(21)*params(24)*params(28)-(1+y(15)*(1-params(23)))*y(8)/y(9);
T(16) = 1/params(1);
T(17) = T(16)*(1+(y(9)-y(1))/y(1))^params(14);
T(18) = (y(22)/params(32))^params(15);
T(19) = T(17)*T(18);
T(20) = ((1+y(15))/(1+params(33)))^params(16);
T(21) = y(17)^(1-params(11));
T(22) = y(16)^params(11);
T(23) = T(21)*T(22);
T(24) = params(11)^params(11);
T(25) = (1-params(11))^(1-params(11));
T(26) = y(9)*exp(y(12))*T(24)*T(25);
T(27) = y(1)/y(9);
T(28) = (-(1+y(14)))/y(14);
T(29) = T(27)^T(28);
T(30) = params(13)*T(29);
T(31) = y(10)/y(9);
T(32) = T(31)^T(28);
T(33) = (1-params(13))*T(32);
T(34) = (-1)/y(14);
T(35) = y(1)^T(34);
T(36) = y(10)^T(34);
T(37) = params(13)*T(35)+(1-params(13))*T(36);
T(38) = T(37)^(-y(14));

end
