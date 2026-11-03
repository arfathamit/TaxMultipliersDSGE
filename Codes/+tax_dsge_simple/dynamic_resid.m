function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
% function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = tax_dsge_simple.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(21, 1);
lhs = y(20)+y(16);
rhs = y(13)*y(9)*y(3)+(1-params(25))*y(18)*y(14)+y(19)+y(21);
residual(1) = lhs - rhs;
lhs = y(23);
rhs = T(3);
residual(2) = lhs - rhs;
lhs = T(6);
rhs = (1-params(25))*y(18)*y(23);
residual(3) = lhs - rhs;
lhs = y(24);
rhs = params(1)*y(32)*y(30)*y(29)+params(1)*y(33)*T(8);
residual(4) = lhs - rhs;
lhs = y(23);
rhs = y(24)*(T(11)-T(10)*params(7)*(T(10)-1))-params(1)*T(12)*T(13);
residual(5) = lhs - rhs;
lhs = y(5)/y(7);
rhs = (1+y(11))*y(23)*params(1)/y(28);
residual(6) = lhs - rhs;
lhs = y(13)*y(23);
rhs = y(24)*(params(9)+params(10)*(y(9)-1));
residual(7) = lhs - rhs;
lhs = y(15);
rhs = y(3)*T(14)+y(16)*T(11);
residual(8) = lhs - rhs;
residual(9) = y(12);
lhs = y(22);
rhs = y(21)*(1-params(30))/params(30);
residual(10) = lhs - rhs;
lhs = y(21);
rhs = params(30)*y(14)*params(25)*y(18);
residual(11) = lhs - rhs;
lhs = y(3);
rhs = (y(13)/T(17))^(1/(params(11)-1));
residual(12) = lhs - rhs;
lhs = y(17);
rhs = T(19);
residual(13) = lhs - rhs;
lhs = y(18);
rhs = (1-params(11))*y(17)/y(14);
residual(14) = lhs - rhs;
lhs = y(27);
rhs = T(22)/T(25);
residual(15) = lhs - rhs;
lhs = y(25);
rhs = T(26)*((y(27)-params(36))/params(36)+y(28)-1+params(12)/(1+params(12))*(y(26)-params(12))/params(12))+params(1)*params(13)*y(34);
residual(16) = lhs - rhs;
lhs = y(19);
rhs = y(17)-y(9)*y(13)*y(3)-y(18)*y(14);
residual(17) = lhs - rhs;
lhs = y(28);
rhs = T(31);
residual(18) = lhs - rhs;
lhs = y(8);
rhs = params(17)*y(1)+x(it_, 2);
residual(19) = lhs - rhs;
lhs = log(y(26));
rhs = params(18)*log(y(6))+(1-params(18))*log(params(12))+x(it_, 3);
residual(20) = lhs - rhs;
lhs = y(10);
rhs = params(6)*y(2)+x(it_, 1);
residual(21) = lhs - rhs;

end
