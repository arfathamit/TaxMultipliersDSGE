function T = static_resid_tt(T, y, x, params)
% function T = static_resid_tt(T, y, x, params)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%
% Output:
%   T         [#temp variables by 1]  double   vector of temporary terms
%

assert(length(T) >= 22);

T(1) = params(3)*y(16)^(params(3)*(1-params(2))-1);
T(2) = (1-y(10))^((1-params(2))*(1-params(3)));
T(3) = exp(y(5))*T(1)*T(2)/(1+params(22));
T(4) = (1-params(3))*y(16)^(params(3)*(1-params(2)));
T(5) = (1-y(10))^((1-params(2))*(1-params(3))-1);
T(6) = exp(y(5))*T(4)*T(5);
T(7) = params(10)/2*(y(3)-1)^2;
T(8) = y(20)/y(1)+y(16)*params(22)+y(9)*params(25)*y(10)/y(1)+y(11)*y(3)*params(24)*y(8)+params(26)*y(15)/y(1)-params(27)*y(13)-params(24)*params(29)*y(12)-y(13)*params(24)*params(28)-(1+y(7)*(1-params(23)))*y(20)/y(1);
T(9) = 1/params(1)*(y(14)/params(32))^params(15);
T(10) = ((1+y(7))/(1+params(33)))^params(16);
T(11) = y(9)^(1-params(11));
T(12) = y(8)^params(11);
T(13) = T(11)*T(12);
T(14) = params(11)^params(11);
T(15) = (1-params(11))^(1-params(11));
T(16) = (y(2)/y(1))^((-(1+y(6)))/y(6));
T(17) = (1-params(13))*T(16);
T(18) = (-1)/y(6);
T(19) = y(1)^T(18);
T(20) = y(2)^T(18);
T(21) = params(13)*T(19)+(1-params(13))*T(20);
T(22) = T(21)^(-y(6));

end
