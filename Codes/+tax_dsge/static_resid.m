function residual = static_resid(T, y, x, params, T_flag)
% function residual = static_resid(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = tax_dsge.static_resid_tt(T, y, x, params);
end
residual = zeros(22, 1);
lhs = y(16);
rhs = ((1+y(7)*(1-params(23)))*y(20)/y(1)+(1-params(25))*y(9)*y(10)/y(1)+(1-params(24))*y(8)*y(3)*y(11)/y(1)+params(24)*params(29)*y(12)+params(27)*y(13)+y(13)*params(24)*params(28)+(1-params(26))*y(15)/y(1)+y(21)-y(13)-y(20)/y(1))/(1+params(22));
residual(1) = lhs - rhs;
lhs = y(17);
rhs = T(3);
residual(2) = lhs - rhs;
lhs = T(6);
rhs = (1-params(25))*y(9)*y(17)/y(1);
residual(3) = lhs - rhs;
lhs = y(18);
rhs = params(1)*(y(3)*y(8)*(1-params(24))*y(17)/y(1)+y(18)*(1-params(8)-params(9)*(y(3)-1)-T(7)));
residual(4) = lhs - rhs;
lhs = y(17)*(1-params(27)-params(24)*params(28));
rhs = y(18)+params(1)*params(29)*params(24)*y(17)*(1-params(28));
residual(5) = lhs - rhs;
lhs = y(17)/y(1);
rhs = (1+y(7)*(1-params(23)))*y(17)*params(1)/y(1);
residual(6) = lhs - rhs;
lhs = y(11)*y(8)*(1-params(24))*y(17)/y(1);
rhs = y(11)*y(18)*(params(9)+(y(3)-1)*params(10));
residual(7) = lhs - rhs;
lhs = y(11);
rhs = y(13)+y(11)*(1-(T(7)+params(8)+params(9)*(y(3)-1)));
residual(8) = lhs - rhs;
lhs = y(12);
rhs = y(12)*(1-params(29))+y(13)*(1-params(28));
residual(9) = lhs - rhs;
lhs = y(14);
rhs = y(16)+y(13)+y(19);
residual(10) = lhs - rhs;
lhs = y(19);
rhs = (1-params(30))*T(8);
residual(11) = lhs - rhs;
lhs = y(21);
rhs = params(30)*T(8);
residual(12) = lhs - rhs;
lhs = y(7);
rhs = (T(9)*T(10)-1)/(1-params(23));
residual(13) = lhs - rhs;
lhs = y(3)*y(11);
rhs = y(1)*y(14)*params(11)/y(8);
residual(14) = lhs - rhs;
lhs = y(10);
rhs = y(1)*y(14)*(1-params(11))/y(9);
residual(15) = lhs - rhs;
lhs = y(22);
rhs = T(13)/(y(1)*exp(y(4))*T(14)*T(15));
residual(16) = lhs - rhs;
lhs = y(2);
rhs = (1+params(12))*params(36)*(1-params(1)*params(13))*((y(22)-params(36))/params(36)+y(1)-1+params(12)/(1+params(12))*(y(6)-params(12))/params(12))+y(2)*params(1)*params(13);
residual(17) = lhs - rhs;
lhs = y(15);
rhs = y(14)*(params(13)*(y(1)-y(22))+T(17)*(y(2)-y(22)));
residual(18) = lhs - rhs;
lhs = y(1);
rhs = T(22);
residual(19) = lhs - rhs;
lhs = y(4);
rhs = y(4)*params(17)+x(3);
residual(20) = lhs - rhs;
lhs = log(y(6));
rhs = log(y(6))*params(18)+(1-params(18))*log(params(12))+x(2);
residual(21) = lhs - rhs;
lhs = y(5);
rhs = y(5)*params(6)+x(1);
residual(22) = lhs - rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
end
