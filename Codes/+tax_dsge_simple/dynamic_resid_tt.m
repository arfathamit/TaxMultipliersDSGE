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

assert(length(T) >= 31);

T(1) = params(3)*y(20)^(params(3)*(1-params(2))-1);
T(2) = (1-y(14))^((1-params(2))*(1-params(3)));
T(3) = exp(y(10))*T(1)*T(2);
T(4) = (1-params(3))*y(20)^(params(3)*(1-params(2)));
T(5) = (1-y(14))^((1-params(2))*(1-params(3))-1);
T(6) = exp(y(10))*T(4)*T(5);
T(7) = params(10)/2;
T(8) = 1-(params(8)+params(9)*(y(29)-1)+T(7)*(y(29)-1)^2);
T(9) = params(7)/2;
T(10) = y(16)/y(4);
T(11) = 1-T(9)*(T(10)-1)^2;
T(12) = y(33)*params(7)*(y(31)/y(16)-1);
T(13) = (y(31)/y(16))^2;
T(14) = 1-(params(8)+params(9)*(y(9)-1)+T(7)*(y(9)-1)^2);
T(15) = params(11)*exp(y(8))*y(9)^params(11);
T(16) = y(14)^(1-params(11));
T(17) = T(15)*T(16);
T(18) = exp(y(8))*(y(9)*y(3))^params(11);
T(19) = T(16)*T(18);
T(20) = y(18)^(1-params(11));
T(21) = y(13)^params(11);
T(22) = T(20)*T(21);
T(23) = params(11)^params(11);
T(24) = (1-params(11))^(1-params(11));
T(25) = y(28)*exp(y(8))*T(23)*T(24);
T(26) = (1+params(12))/params(12)*params(36)*(1-params(1)*params(13));
T(27) = (-1)/y(26);
T(28) = y(7)^T(27);
T(29) = exp(y(25))^T(27);
T(30) = params(13)*T(28)+(1-params(13))*T(29);
T(31) = T(30)^(-y(26));

end
