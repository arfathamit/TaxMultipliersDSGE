function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
% function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
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
%   g1
%

if T_flag
    T = tax_dsge_simple.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(21, 37);
g1(1,9)=(-(y(13)*y(3)));
g1(1,13)=(-(y(9)*y(3)));
g1(1,14)=(-((1-params(25))*y(18)));
g1(1,3)=(-(y(13)*y(9)));
g1(1,16)=1;
g1(1,18)=(-((1-params(25))*y(14)));
g1(1,19)=(-1);
g1(1,20)=1;
g1(1,21)=(-1);
g1(2,10)=(-T(3));
g1(2,14)=T(43);
g1(2,20)=T(60);
g1(2,23)=1;
g1(3,10)=T(6);
g1(3,14)=T(45);
g1(3,18)=(-((1-params(25))*y(23)));
g1(3,20)=T(62);
g1(3,23)=(-((1-params(25))*y(18)));
g1(4,29)=(-(params(1)*y(32)*y(30)+params(1)*y(33)*(-(params(9)+T(7)*2*(y(29)-1)))));
g1(4,30)=(-(params(1)*y(32)*y(29)));
g1(4,32)=(-(y(29)*params(1)*y(30)));
g1(4,24)=1;
g1(4,33)=(-(params(1)*T(8)));
g1(5,4)=(-(y(24)*(T(49)-(params(7)*(T(10)-1)*T(48)+T(10)*params(7)*T(48)))));
g1(5,16)=(-(y(24)*(T(51)-(params(7)*(T(10)-1)*T(50)+T(10)*params(7)*T(50)))-params(1)*(T(13)*y(33)*params(7)*T(52)+T(12)*T(54))));
g1(5,31)=params(1)*(T(13)*y(33)*params(7)*T(55)+T(12)*T(56));
g1(5,23)=1;
g1(5,24)=(-(T(11)-T(10)*params(7)*(T(10)-1)));
g1(5,33)=params(1)*T(13)*params(7)*(y(31)/y(16)-1);
g1(6,11)=(-(y(23)*params(1)/y(28)));
g1(6,5)=1/y(7);
g1(6,23)=(-((1+y(11))*params(1)/y(28)));
g1(6,7)=(-y(5))/(y(7)*y(7));
g1(6,28)=(-((-((1+y(11))*y(23)*params(1)))/(y(28)*y(28))));
g1(7,9)=(-(y(24)*params(10)));
g1(7,13)=y(23);
g1(7,23)=y(13);
g1(7,24)=(-(params(9)+params(10)*(y(9)-1)));
g1(8,9)=(-(y(3)*(-(params(9)+T(7)*2*(y(9)-1)))));
g1(8,3)=(-T(14));
g1(8,15)=1;
g1(8,4)=(-(y(16)*T(49)));
g1(8,16)=(-(T(11)+y(16)*T(51)));
g1(9,12)=1;
g1(10,21)=(-((1-params(30))/params(30)));
g1(10,22)=1;
g1(11,14)=(-(params(30)*params(25)*y(18)));
g1(11,18)=(-(params(30)*params(25)*y(14)));
g1(11,21)=1;
g1(12,8)=(-(T(32)*T(33)));
g1(12,9)=(-(T(33)*T(37)));
g1(12,13)=(-(T(33)*T(39)));
g1(12,14)=(-(T(33)*T(47)));
g1(12,3)=1;
g1(13,8)=(-T(19));
g1(13,9)=(-(T(16)*exp(y(8))*y(3)*T(38)));
g1(13,14)=(-(T(18)*T(46)));
g1(13,3)=(-(T(16)*exp(y(8))*y(9)*T(38)));
g1(13,17)=1;
g1(14,14)=(-((1-params(11))*(-y(17))/(y(14)*y(14))));
g1(14,17)=(-((1-params(11))*1/y(14)));
g1(14,18)=1;
g1(15,8)=(-((-(T(22)*T(25)))/(T(25)*T(25))));
g1(15,13)=(-(T(41)/T(25)));
g1(15,18)=(-(T(58)/T(25)));
g1(15,27)=1;
g1(15,28)=(-((-(T(22)*T(24)*exp(y(8))*T(23)))/(T(25)*T(25))));
g1(16,25)=1;
g1(16,34)=(-(params(1)*params(13)));
g1(16,26)=(-(T(26)*params(12)/(1+params(12))*1/params(12)));
g1(16,27)=(-(T(26)*1/params(36)));
g1(16,28)=(-T(26));
g1(17,9)=y(13)*y(3);
g1(17,13)=y(9)*y(3);
g1(17,14)=y(18);
g1(17,3)=y(13)*y(9);
g1(17,17)=(-1);
g1(17,18)=y(14);
g1(17,19)=1;
g1(18,25)=(-((1-params(13))*T(63)*T(64)));
g1(18,26)=(-(T(31)*T(69)));
g1(18,7)=(-(T(64)*T(71)));
g1(18,28)=1;
g1(19,1)=(-params(17));
g1(19,8)=1;
g1(19,36)=(-1);
g1(20,6)=(-(params(18)*1/y(6)));
g1(20,26)=1/y(26);
g1(20,37)=(-1);
g1(21,2)=(-params(6));
g1(21,10)=1;
g1(21,35)=(-1);

end
