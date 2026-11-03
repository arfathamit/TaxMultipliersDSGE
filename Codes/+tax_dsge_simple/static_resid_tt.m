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

assert(length(T) >= 24);

T(1) = params(3)*y(13)^(params(3)*(1-params(2))-1);
T(2) = (1-y(7))^((1-params(2))*(1-params(3)));
T(3) = exp(y(3))*T(1)*T(2);
T(4) = (1-params(3))*y(13)^(params(3)*(1-params(2)));
T(5) = (1-y(7))^((1-params(2))*(1-params(3))-1);
T(6) = exp(y(3))*T(4)*T(5);
T(7) = 1-(params(8)+params(9)*(y(2)-1)+params(10)/2*(y(2)-1)^2);
T(8) = params(11)*exp(y(1))*y(2)^params(11);
T(9) = y(7)^(1-params(11));
T(10) = T(8)*T(9);
T(11) = exp(y(1))*(y(2)*y(8))^params(11);
T(12) = T(9)*T(11);
T(13) = y(11)^(1-params(11));
T(14) = y(6)^params(11);
T(15) = T(13)*T(14);
T(16) = params(11)^params(11);
T(17) = (1-params(11))^(1-params(11));
T(18) = y(21)*exp(y(1))*T(16)*T(17);
T(19) = (1+params(12))/params(12)*params(36)*(1-params(1)*params(13));
T(20) = (-1)/y(19);
T(21) = y(21)^T(20);
T(22) = exp(y(18))^T(20);
T(23) = params(13)*T(21)+(1-params(13))*T(22);
T(24) = T(23)^(-y(19));

end
