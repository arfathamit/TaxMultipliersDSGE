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
    T = tax_dsge_simple.static_resid_tt(T, y, x, params);
end
residual = zeros(21, 1);
lhs = y(13)+y(9);
rhs = y(6)*y(2)*y(8)+(1-params(25))*y(11)*y(7)+y(12)+y(14);
residual(1) = lhs - rhs;
lhs = y(16);
rhs = T(3);
residual(2) = lhs - rhs;
lhs = T(6);
rhs = (1-params(25))*y(11)*y(16);
residual(3) = lhs - rhs;
lhs = y(17);
rhs = y(2)*y(6)*y(16)*params(1)+y(17)*params(1)*T(7);
residual(4) = lhs - rhs;
lhs = y(16);
rhs = y(17);
residual(5) = lhs - rhs;
lhs = y(16)/y(21);
rhs = y(16)*params(1)*(1+y(4))/y(21);
residual(6) = lhs - rhs;
lhs = y(6)*y(16);
rhs = y(17)*(params(9)+(y(2)-1)*params(10));
residual(7) = lhs - rhs;
lhs = y(8);
rhs = y(9)+y(8)*T(7);
residual(8) = lhs - rhs;
residual(9) = y(5);
lhs = y(15);
rhs = y(14)*(1-params(30))/params(30);
residual(10) = lhs - rhs;
lhs = y(14);
rhs = params(30)*y(7)*params(25)*y(11);
residual(11) = lhs - rhs;
lhs = y(8);
rhs = (y(6)/T(10))^(1/(params(11)-1));
residual(12) = lhs - rhs;
lhs = y(10);
rhs = T(12);
residual(13) = lhs - rhs;
lhs = y(11);
rhs = (1-params(11))*y(10)/y(7);
residual(14) = lhs - rhs;
lhs = y(20);
rhs = T(15)/T(18);
residual(15) = lhs - rhs;
lhs = y(18);
rhs = T(19)*((y(20)-params(36))/params(36)+y(21)-1+params(12)/(1+params(12))*(y(19)-params(12))/params(12))+y(18)*params(1)*params(13);
residual(16) = lhs - rhs;
lhs = y(12);
rhs = y(10)-y(2)*y(6)*y(8)-y(11)*y(7);
residual(17) = lhs - rhs;
lhs = y(21);
rhs = T(24);
residual(18) = lhs - rhs;
lhs = y(1);
rhs = y(1)*params(17)+x(2);
residual(19) = lhs - rhs;
lhs = log(y(19));
rhs = log(y(19))*params(18)+(1-params(18))*log(params(12))+x(3);
residual(20) = lhs - rhs;
lhs = y(3);
rhs = y(3)*params(6)+x(1);
residual(21) = lhs - rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
end
