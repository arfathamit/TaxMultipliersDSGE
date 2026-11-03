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
    T = tax_dsge.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(22, 41);
g1(1,9)=(-(((-((1+y(15)*(1-params(23)))*y(8)))/(y(9)*y(9))+(-((1-params(25))*y(17)*y(18)))/(y(9)*y(9))+(-((1-params(24))*y(16)*y(11)*y(5)))/(y(9)*y(9))+(-((1-params(26))*y(23)))/(y(9)*y(9))-(-y(28))/(y(9)*y(9)))/(1+params(22))));
g1(1,11)=(-((1-params(24))*y(16)*y(5)/y(9)/(1+params(22))));
g1(1,15)=(-((1-params(23))*y(8)/y(9)/(1+params(22))));
g1(1,16)=(-(y(5)*(1-params(24))*y(11)/y(9)/(1+params(22))));
g1(1,17)=(-((1-params(25))*y(18)/y(9)/(1+params(22))));
g1(1,18)=(-((1-params(25))*y(17)/y(9)/(1+params(22))));
g1(1,5)=(-((1-params(24))*y(16)*y(11)/y(9)/(1+params(22))));
g1(1,6)=(-(params(24)*params(29)/(1+params(22))));
g1(1,21)=(-((params(27)+params(24)*params(28)-1)/(1+params(22))));
g1(1,23)=(-((1-params(26))/y(9)/(1+params(22))));
g1(1,24)=1;
g1(1,8)=(-((1+y(15)*(1-params(23)))/y(9)/(1+params(22))));
g1(1,28)=(-((-T(41))/(1+params(22))));
g1(1,29)=(-(1/(1+params(22))));
g1(2,13)=(-T(3));
g1(2,18)=T(60);
g1(2,24)=T(74);
g1(2,25)=1;
g1(3,9)=(-((-((1-params(25))*y(17)*y(25)))/(y(9)*y(9))));
g1(3,13)=T(6);
g1(3,17)=(-((1-params(25))*y(25)/y(9)));
g1(3,18)=T(62);
g1(3,24)=T(76);
g1(3,25)=(-((1-params(25))*y(17)/y(9)));
g1(4,31)=(-(params(1)*(-((1-params(24))*y(37)*y(35)*y(33)))/(y(31)*y(31))));
g1(4,33)=(-(params(1)*((1-params(24))*y(37)*y(35)/y(31)+y(38)*((-params(9))-T(7)*2*(y(33)-1)))));
g1(4,35)=(-(params(1)*(1-params(24))*y(37)*y(33)/y(31)));
g1(4,37)=(-(params(1)*y(33)*(1-params(24))*y(35)/y(31)));
g1(4,26)=1;
g1(4,38)=(-(params(1)*T(8)));
g1(5,7)=(-(y(26)*(T(64)-(params(7)*(T(10)-1)*T(63)+T(10)*params(7)*T(63)))));
g1(5,21)=(-(y(26)*(T(66)-(params(7)*(T(10)-1)*T(65)+T(10)*params(7)*T(65)))+params(1)*(T(13)*y(38)*params(7)*T(67)+T(12)*T(69))));
g1(5,36)=(-(params(1)*(T(13)*y(38)*params(7)*T(70)+T(12)*T(71))));
g1(5,25)=1-params(27)-params(24)*params(28);
g1(5,37)=(-(params(1)*params(24)*params(29)*(1-params(28))));
g1(5,26)=(-(T(11)-T(10)*params(7)*(T(10)-1)));
g1(5,38)=(-(params(1)*T(13)*params(7)*(y(36)/y(21)-1)));
g1(6,9)=(-y(25))/(y(9)*y(9));
g1(6,31)=(-((-(params(1)*y(37)*(1+(1-params(23))*y(34))))/(y(31)*y(31))));
g1(6,34)=(-((1-params(23))*params(1)*y(37)/y(31)));
g1(6,25)=T(41);
g1(6,37)=(-(params(1)*(1+(1-params(23))*y(34))/y(31)));
g1(7,9)=(-(y(5)*y(16)*(1-params(24))*y(25)))/(y(9)*y(9));
g1(7,11)=(-(y(5)*y(26)*params(10)));
g1(7,16)=y(5)*(1-params(24))*y(25)/y(9);
g1(7,5)=y(16)*(1-params(24))*y(25)/y(9)-y(26)*(params(9)+params(10)*(y(11)-1));
g1(7,25)=(1-params(24))*y(16)*y(5)/y(9);
g1(7,26)=(-(y(5)*(params(9)+params(10)*(y(11)-1))));
g1(8,11)=(-(y(5)*(-(params(9)+T(7)*2*(y(11)-1)))));
g1(8,5)=(-T(14));
g1(8,19)=1;
g1(8,7)=(-(y(21)*T(64)));
g1(8,21)=(-(T(11)+y(21)*T(66)));
g1(9,6)=(-(1-params(29)));
g1(9,20)=1;
g1(9,21)=(-(1-params(28)));
g1(10,21)=(-1);
g1(10,22)=1;
g1(10,24)=(-1);
g1(10,27)=(-1);
g1(11,9)=(-((1-params(30))*T(45)));
g1(11,11)=(-((1-params(30))*y(5)*params(24)*y(16)));
g1(11,15)=(-((1-params(30))*(-((1-params(23))*y(8)/y(9)))));
g1(11,16)=(-((1-params(30))*y(5)*params(24)*y(11)));
g1(11,17)=(-((1-params(30))*params(25)*y(18)/y(9)));
g1(11,18)=(-((1-params(30))*params(25)*y(17)/y(9)));
g1(11,5)=(-((1-params(30))*y(11)*params(24)*y(16)));
g1(11,6)=(-((1-params(30))*(-(params(24)*params(29)))));
g1(11,21)=(-((1-params(30))*((-params(27))-params(24)*params(28))));
g1(11,23)=(-((1-params(30))*params(26)/y(9)));
g1(11,24)=(-(params(22)*(1-params(30))));
g1(11,27)=1;
g1(11,8)=(-((1-params(30))*(-((1+y(15)*(1-params(23)))/y(9)))));
g1(11,28)=(-((1-params(30))*T(41)));
g1(12,9)=(-(params(30)*T(45)));
g1(12,11)=(-(params(30)*y(5)*params(24)*y(16)));
g1(12,15)=(-(params(30)*(-((1-params(23))*y(8)/y(9)))));
g1(12,16)=(-(params(30)*y(5)*params(24)*y(11)));
g1(12,17)=(-(params(30)*params(25)*y(18)/y(9)));
g1(12,18)=(-(params(30)*params(25)*y(17)/y(9)));
g1(12,5)=(-(params(30)*y(11)*params(24)*y(16)));
g1(12,6)=(-(params(30)*(-(params(24)*params(29)))));
g1(12,21)=(-(params(30)*((-params(27))-params(24)*params(28))));
g1(12,23)=(-(params(30)*params(26)/y(9)));
g1(12,24)=(-(params(22)*params(30)));
g1(12,8)=(-(params(30)*(-((1+y(15)*(1-params(23)))/y(9)))));
g1(12,28)=(-(params(30)*T(41)));
g1(12,29)=1;
g1(13,1)=(-(T(20)*T(18)*T(16)*T(39)*T(40)/(1-params(23))));
g1(13,9)=(-(T(20)*T(18)*T(16)*T(40)*1/y(1)/(1-params(23))));
g1(13,15)=(-(T(19)*T(54)/(1-params(23))));
g1(13,34)=1;
g1(13,22)=(-(T(20)*T(17)*T(72)/(1-params(23))));
g1(14,9)=(-(y(22)*params(11)/y(16)));
g1(14,11)=y(5);
g1(14,16)=(-(y(9)*y(22)*(-params(11))/(y(16)*y(16))));
g1(14,5)=y(11);
g1(14,22)=(-(y(9)*params(11)/y(16)));
g1(15,9)=(-(y(22)*(1-params(11))/y(17)));
g1(15,17)=(-(y(9)*y(22)*(-(1-params(11)))/(y(17)*y(17))));
g1(15,18)=1;
g1(15,22)=(-(y(9)*(1-params(11))/y(17)));
g1(16,9)=(-((-(T(23)*T(25)*exp(y(12))*T(24)))/T(46)));
g1(16,12)=(-((-(T(23)*T(26)))/T(46)));
g1(16,16)=(-(T(56)/T(26)));
g1(16,17)=(-(T(58)/T(26)));
g1(16,30)=1;
g1(17,9)=(-((1+params(12))*params(36)*(1-params(1)*params(13))));
g1(17,10)=1;
g1(17,32)=(-(params(1)*params(13)));
g1(17,14)=(-((1+params(12))*params(36)*(1-params(1)*params(13))*params(12)/(1+params(12))*1/params(12)));
g1(17,30)=(-((1+params(12))*params(36)*(1-params(1)*params(13))*1/params(36)));
g1(18,1)=(-(y(22)*(T(30)+(y(1)-y(30))*params(13)*T(41)*T(42))));
g1(18,9)=(-(y(22)*((y(1)-y(30))*params(13)*T(42)*T(47)+(y(10)-y(30))*(1-params(13))*T(48)*T(49))));
g1(18,10)=(-(y(22)*(T(33)+(y(10)-y(30))*(1-params(13))*T(41)*T(49))));
g1(18,14)=(-(y(22)*((y(1)-y(30))*params(13)*T(29)*T(51)*log(T(27))+(y(10)-y(30))*(1-params(13))*T(32)*T(51)*log(T(31)))));
g1(18,22)=(-(T(30)*(y(1)-y(30))+T(33)*(y(10)-y(30))));
g1(18,23)=1;
g1(18,30)=(-(y(22)*((-T(30))-T(33))));
g1(19,1)=(-(T(43)*T(44)));
g1(19,9)=1;
g1(19,10)=(-(T(44)*T(50)));
g1(19,14)=(-(T(38)*T(53)));
g1(20,2)=(-params(17));
g1(20,12)=1;
g1(20,41)=(-1);
g1(21,4)=(-(params(18)*1/y(4)));
g1(21,14)=1/y(14);
g1(21,40)=(-1);
g1(22,3)=(-params(6));
g1(22,13)=1;
g1(22,39)=(-1);

end
