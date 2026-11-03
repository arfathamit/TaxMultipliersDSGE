function g1 = static_g1(T, y, x, params, T_flag)
% function g1 = static_g1(T, y, x, params, T_flag)
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
%   g1
%

if T_flag
    T = tax_dsge.static_g1_tt(T, y, x, params);
end
g1 = zeros(22, 22);
g1(1,1)=(-(((-((1+y(7)*(1-params(23)))*y(20)))/(y(1)*y(1))+(-((1-params(25))*y(9)*y(10)))/(y(1)*y(1))+(-((1-params(24))*y(8)*y(3)*y(11)))/(y(1)*y(1))+(-((1-params(26))*y(15)))/(y(1)*y(1))-(-y(20))/(y(1)*y(1)))/(1+params(22))));
g1(1,3)=(-((1-params(24))*y(8)*y(11)/y(1)/(1+params(22))));
g1(1,7)=(-((1-params(23))*y(20)/y(1)/(1+params(22))));
g1(1,8)=(-(y(11)*(1-params(24))*y(3)/y(1)/(1+params(22))));
g1(1,9)=(-((1-params(25))*y(10)/y(1)/(1+params(22))));
g1(1,10)=(-((1-params(25))*y(9)/y(1)/(1+params(22))));
g1(1,11)=(-((1-params(24))*y(8)*y(3)/y(1)/(1+params(22))));
g1(1,12)=(-(params(24)*params(29)/(1+params(22))));
g1(1,13)=(-((params(27)+params(24)*params(28)-1)/(1+params(22))));
g1(1,15)=(-((1-params(26))/y(1)/(1+params(22))));
g1(1,16)=1;
g1(1,20)=(-(((1+y(7)*(1-params(23)))/y(1)-T(26))/(1+params(22))));
g1(1,21)=(-(1/(1+params(22))));
g1(2,5)=(-T(3));
g1(2,10)=(-(exp(y(5))*T(1)*(-(getPowerDeriv(1-y(10),(1-params(2))*(1-params(3)),1)))/(1+params(22))));
g1(2,16)=(-(exp(y(5))*T(2)*params(3)*getPowerDeriv(y(16),params(3)*(1-params(2))-1,1)/(1+params(22))));
g1(2,17)=1;
g1(3,1)=(-((-((1-params(25))*y(9)*y(17)))/(y(1)*y(1))));
g1(3,5)=T(6);
g1(3,9)=(-((1-params(25))*y(17)/y(1)));
g1(3,10)=exp(y(5))*T(4)*(-(getPowerDeriv(1-y(10),(1-params(2))*(1-params(3))-1,1)));
g1(3,16)=exp(y(5))*T(5)*(1-params(3))*getPowerDeriv(y(16),params(3)*(1-params(2)),1);
g1(3,17)=(-((1-params(25))*y(9)/y(1)));
g1(4,1)=(-(params(1)*(-(y(3)*y(8)*(1-params(24))*y(17)))/(y(1)*y(1))));
g1(4,3)=(-(params(1)*(y(8)*(1-params(24))*y(17)/y(1)+y(18)*((-params(9))-params(10)/2*2*(y(3)-1)))));
g1(4,8)=(-(params(1)*y(3)*(1-params(24))*y(17)/y(1)));
g1(4,17)=(-(params(1)*(1-params(24))*y(8)*y(3)/y(1)));
g1(4,18)=1-params(1)*(1-params(8)-params(9)*(y(3)-1)-T(7));
g1(5,17)=1-params(27)-params(24)*params(28)-params(1)*params(24)*params(29)*(1-params(28));
g1(5,18)=(-1);
g1(6,1)=(-y(17))/(y(1)*y(1))-(-((1+y(7)*(1-params(23)))*y(17)*params(1)))/(y(1)*y(1));
g1(6,7)=(-((1-params(23))*y(17)*params(1)/y(1)));
g1(6,17)=T(26)-(1+y(7)*(1-params(23)))*params(1)/y(1);
g1(7,1)=(-(y(11)*y(8)*(1-params(24))*y(17)))/(y(1)*y(1));
g1(7,3)=(-(y(11)*y(18)*params(10)));
g1(7,8)=y(11)*(1-params(24))*y(17)/y(1);
g1(7,11)=y(8)*(1-params(24))*y(17)/y(1)-y(18)*(params(9)+(y(3)-1)*params(10));
g1(7,17)=(1-params(24))*y(8)*y(11)/y(1);
g1(7,18)=(-(y(11)*(params(9)+(y(3)-1)*params(10))));
g1(8,3)=(-(y(11)*(-(params(9)+params(10)/2*2*(y(3)-1)))));
g1(8,11)=1-(1-(T(7)+params(8)+params(9)*(y(3)-1)));
g1(8,13)=(-1);
g1(9,12)=1-(1-params(29));
g1(9,13)=(-(1-params(28)));
g1(10,13)=(-1);
g1(10,14)=1;
g1(10,16)=(-1);
g1(10,19)=(-1);
g1(11,1)=(-((1-params(30))*T(23)));
g1(11,3)=(-((1-params(30))*y(11)*params(24)*y(8)));
g1(11,7)=(-((1-params(30))*(-((1-params(23))*y(20)/y(1)))));
g1(11,8)=(-((1-params(30))*y(11)*params(24)*y(3)));
g1(11,9)=(-((1-params(30))*params(25)*y(10)/y(1)));
g1(11,10)=(-((1-params(30))*params(25)*y(9)/y(1)));
g1(11,11)=(-((1-params(30))*y(3)*params(24)*y(8)));
g1(11,12)=(-((1-params(30))*(-(params(24)*params(29)))));
g1(11,13)=(-((1-params(30))*((-params(27))-params(24)*params(28))));
g1(11,15)=(-((1-params(30))*params(26)/y(1)));
g1(11,16)=(-(params(22)*(1-params(30))));
g1(11,19)=1;
g1(11,20)=(-((1-params(30))*T(27)));
g1(12,1)=(-(params(30)*T(23)));
g1(12,3)=(-(params(30)*y(11)*params(24)*y(8)));
g1(12,7)=(-(params(30)*(-((1-params(23))*y(20)/y(1)))));
g1(12,8)=(-(params(30)*y(11)*params(24)*y(3)));
g1(12,9)=(-(params(30)*params(25)*y(10)/y(1)));
g1(12,10)=(-(params(30)*params(25)*y(9)/y(1)));
g1(12,11)=(-(params(30)*y(3)*params(24)*y(8)));
g1(12,12)=(-(params(30)*(-(params(24)*params(29)))));
g1(12,13)=(-(params(30)*((-params(27))-params(24)*params(28))));
g1(12,15)=(-(params(30)*params(26)/y(1)));
g1(12,16)=(-(params(22)*params(30)));
g1(12,20)=(-(params(30)*T(27)));
g1(12,21)=1;
g1(13,7)=1-T(9)*1/(1+params(33))*getPowerDeriv((1+y(7))/(1+params(33)),params(16),1)/(1-params(23));
g1(13,14)=(-(T(10)*1/params(1)*1/params(32)*getPowerDeriv(y(14)/params(32),params(15),1)/(1-params(23))));
g1(14,1)=(-(y(14)*params(11)/y(8)));
g1(14,3)=y(11);
g1(14,8)=(-(y(1)*y(14)*(-params(11))/(y(8)*y(8))));
g1(14,11)=y(3);
g1(14,14)=(-(y(1)*params(11)/y(8)));
g1(15,1)=(-(y(14)*(1-params(11))/y(9)));
g1(15,9)=(-(y(1)*y(14)*(-(1-params(11)))/(y(9)*y(9))));
g1(15,10)=1;
g1(15,14)=(-(y(1)*(1-params(11))/y(9)));
g1(16,1)=(-((-(T(13)*T(15)*exp(y(4))*T(14)))/(y(1)*exp(y(4))*T(14)*T(15)*y(1)*exp(y(4))*T(14)*T(15))));
g1(16,4)=(-((-(T(13)*y(1)*exp(y(4))*T(14)*T(15)))/(y(1)*exp(y(4))*T(14)*T(15)*y(1)*exp(y(4))*T(14)*T(15))));
g1(16,8)=(-(T(11)*getPowerDeriv(y(8),params(11),1)/(y(1)*exp(y(4))*T(14)*T(15))));
g1(16,9)=(-(T(12)*getPowerDeriv(y(9),1-params(11),1)/(y(1)*exp(y(4))*T(14)*T(15))));
g1(16,22)=1;
g1(17,1)=(-((1+params(12))*params(36)*(1-params(1)*params(13))));
g1(17,2)=1-params(1)*params(13);
g1(17,6)=(-((1+params(12))*params(36)*(1-params(1)*params(13))*params(12)/(1+params(12))*1/params(12)));
g1(17,22)=(-((1+params(12))*params(36)*(1-params(1)*params(13))*1/params(36)));
g1(18,1)=(-(y(14)*(params(13)+(y(2)-y(22))*(1-params(13))*(-y(2))/(y(1)*y(1))*T(24))));
g1(18,2)=(-(y(14)*(T(17)+(y(2)-y(22))*(1-params(13))*T(24)*T(26))));
g1(18,6)=(-(y(14)*(y(2)-y(22))*(1-params(13))*T(16)*1/(y(6)*y(6))*log(y(2)/y(1))));
g1(18,14)=(-(params(13)*(y(1)-y(22))+T(17)*(y(2)-y(22))));
g1(18,15)=1;
g1(18,22)=(-(y(14)*((-params(13))-T(17))));
g1(19,1)=1-params(13)*getPowerDeriv(y(1),T(18),1)*T(25);
g1(19,2)=(-(T(25)*(1-params(13))*getPowerDeriv(y(2),T(18),1)));
g1(19,6)=(-(T(22)*((-y(6))*(params(13)*T(19)*1/(y(6)*y(6))*log(y(1))+(1-params(13))*T(20)*1/(y(6)*y(6))*log(y(2)))/T(21)-log(T(21)))));
g1(20,4)=1-params(17);
g1(21,6)=1/y(6)-params(18)*1/y(6);
g1(22,5)=1-params(6);
if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
end
end
