function T = static_g1_tt(T, y, x, params)
% function T = static_g1_tt(T, y, x, params)
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

assert(length(T) >= 27);

T = tax_dsge.static_resid_tt(T, y, x, params);

T(23) = (-y(20))/(y(1)*y(1))+(-(y(9)*params(25)*y(10)))/(y(1)*y(1))+(-(params(26)*y(15)))/(y(1)*y(1))-(-((1+y(7)*(1-params(23)))*y(20)))/(y(1)*y(1));
T(24) = getPowerDeriv(y(2)/y(1),(-(1+y(6)))/y(6),1);
T(25) = getPowerDeriv(T(21),(-y(6)),1);
T(26) = 1/y(1);
T(27) = T(26)-(1+y(7)*(1-params(23)))/y(1);

end
