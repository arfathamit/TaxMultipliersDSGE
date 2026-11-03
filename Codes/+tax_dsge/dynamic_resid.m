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
    T = tax_dsge.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(22, 1);
lhs = y(24);
rhs = ((1+y(15)*(1-params(23)))*y(8)/y(9)+(1-params(25))*y(17)*y(18)/y(9)+(1-params(24))*y(16)*y(11)*y(5)/y(9)+params(24)*params(29)*y(6)+params(27)*y(21)+y(21)*params(24)*params(28)+(1-params(26))*y(23)/y(9)+y(29)-y(21)-y(28)/y(9))/(1+params(22));
residual(1) = lhs - rhs;
lhs = y(25);
rhs = T(3);
residual(2) = lhs - rhs;
lhs = T(6);
rhs = (1-params(25))*y(17)*y(25)/y(9);
residual(3) = lhs - rhs;
lhs = y(26);
rhs = params(1)*((1-params(24))*y(37)*y(35)*y(33)/y(31)+y(38)*T(8));
residual(4) = lhs - rhs;
lhs = y(25)*(1-params(27)-params(24)*params(28));
rhs = y(26)*(T(11)-T(10)*params(7)*(T(10)-1))+params(1)*(params(29)*params(24)*y(37)*(1-params(28))+T(12)*T(13));
residual(5) = lhs - rhs;
lhs = y(25)/y(9);
rhs = params(1)*y(37)*(1+(1-params(23))*y(34))/y(31);
residual(6) = lhs - rhs;
lhs = y(5)*y(16)*(1-params(24))*y(25)/y(9);
rhs = y(5)*y(26)*(params(9)+params(10)*(y(11)-1));
residual(7) = lhs - rhs;
lhs = y(19);
rhs = y(5)*T(14)+y(21)*T(11);
residual(8) = lhs - rhs;
lhs = y(20);
rhs = y(6)*(1-params(29))+y(21)*(1-params(28));
residual(9) = lhs - rhs;
lhs = y(22);
rhs = y(24)+y(21)+y(27);
residual(10) = lhs - rhs;
lhs = y(27);
rhs = (1-params(30))*T(15);
residual(11) = lhs - rhs;
lhs = y(29);
rhs = params(30)*T(15);
residual(12) = lhs - rhs;
lhs = y(34);
rhs = (T(19)*T(20)-1)/(1-params(23));
residual(13) = lhs - rhs;
lhs = y(11)*y(5);
rhs = y(9)*y(22)*params(11)/y(16);
residual(14) = lhs - rhs;
lhs = y(18);
rhs = y(9)*y(22)*(1-params(11))/y(17);
residual(15) = lhs - rhs;
lhs = y(30);
rhs = T(23)/T(26);
residual(16) = lhs - rhs;
lhs = y(10);
rhs = (1+params(12))*params(36)*(1-params(1)*params(13))*((y(30)-params(36))/params(36)+y(9)-1+params(12)/(1+params(12))*(y(14)-params(12))/params(12))+params(1)*params(13)*y(32);
residual(17) = lhs - rhs;
lhs = y(23);
rhs = y(22)*(T(30)*(y(1)-y(30))+T(33)*(y(10)-y(30)));
residual(18) = lhs - rhs;
lhs = y(9);
rhs = T(38);
residual(19) = lhs - rhs;
lhs = y(12);
rhs = params(17)*y(2)+x(it_, 3);
residual(20) = lhs - rhs;
lhs = log(y(14));
rhs = params(18)*log(y(4))+(1-params(18))*log(params(12))+x(it_, 2);
residual(21) = lhs - rhs;
lhs = y(13);
rhs = params(6)*y(3)+x(it_, 1);
residual(22) = lhs - rhs;

end
