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
    T = tax_dsge_simple.static_g1_tt(T, y, x, params);
end
g1 = zeros(21, 21);
g1(1,2)=(-(y(6)*y(8)));
g1(1,6)=(-(y(2)*y(8)));
g1(1,7)=(-((1-params(25))*y(11)));
g1(1,8)=(-(y(6)*y(2)));
g1(1,9)=1;
g1(1,11)=(-((1-params(25))*y(7)));
g1(1,12)=(-1);
g1(1,13)=1;
g1(1,14)=(-1);
g1(2,3)=(-T(3));
g1(2,7)=(-(exp(y(3))*T(1)*(-(getPowerDeriv(1-y(7),(1-params(2))*(1-params(3)),1)))));
g1(2,13)=(-(exp(y(3))*T(2)*params(3)*getPowerDeriv(y(13),params(3)*(1-params(2))-1,1)));
g1(2,16)=1;
g1(3,3)=T(6);
g1(3,7)=exp(y(3))*T(4)*(-(getPowerDeriv(1-y(7),(1-params(2))*(1-params(3))-1,1)));
g1(3,11)=(-((1-params(25))*y(16)));
g1(3,13)=exp(y(3))*T(5)*(1-params(3))*getPowerDeriv(y(13),params(3)*(1-params(2)),1);
g1(3,16)=(-((1-params(25))*y(11)));
g1(4,2)=(-(y(6)*y(16)*params(1)+y(17)*params(1)*(-(params(9)+params(10)/2*2*(y(2)-1)))));
g1(4,6)=(-(y(2)*y(16)*params(1)));
g1(4,16)=(-(y(2)*y(6)*params(1)));
g1(4,17)=1-params(1)*T(7);
g1(5,16)=1;
g1(5,17)=(-1);
g1(6,4)=(-(y(16)*params(1)/y(21)));
g1(6,16)=1/y(21)-params(1)*(1+y(4))/y(21);
g1(6,21)=(-y(16))/(y(21)*y(21))-(-(y(16)*params(1)*(1+y(4))))/(y(21)*y(21));
g1(7,2)=(-(y(17)*params(10)));
g1(7,6)=y(16);
g1(7,16)=y(6);
g1(7,17)=(-(params(9)+(y(2)-1)*params(10)));
g1(8,2)=(-(y(8)*(-(params(9)+params(10)/2*2*(y(2)-1)))));
g1(8,8)=1-T(7);
g1(8,9)=(-1);
g1(9,5)=1;
g1(10,14)=(-((1-params(30))/params(30)));
g1(10,15)=1;
g1(11,7)=(-(params(30)*params(25)*y(11)));
g1(11,11)=(-(params(30)*params(25)*y(7)));
g1(11,14)=1;
g1(12,1)=(-((-(y(6)*T(10)))/(T(10)*T(10))*T(25)));
g1(12,2)=(-(T(25)*(-(y(6)*T(9)*params(11)*exp(y(1))*getPowerDeriv(y(2),params(11),1)))/(T(10)*T(10))));
g1(12,6)=(-(T(25)*1/T(10)));
g1(12,7)=(-(T(25)*(-(y(6)*T(8)*T(27)))/(T(10)*T(10))));
g1(12,8)=1;
g1(13,1)=(-T(12));
g1(13,2)=(-(T(9)*exp(y(1))*y(8)*T(26)));
g1(13,7)=(-(T(11)*T(27)));
g1(13,8)=(-(T(9)*exp(y(1))*y(2)*T(26)));
g1(13,10)=1;
g1(14,7)=(-((1-params(11))*(-y(10))/(y(7)*y(7))));
g1(14,10)=(-((1-params(11))*1/y(7)));
g1(14,11)=1;
g1(15,1)=(-((-(T(15)*T(18)))/(T(18)*T(18))));
g1(15,6)=(-(T(13)*getPowerDeriv(y(6),params(11),1)/T(18)));
g1(15,11)=(-(T(14)*getPowerDeriv(y(11),1-params(11),1)/T(18)));
g1(15,20)=1;
g1(15,21)=(-((-(T(15)*T(17)*exp(y(1))*T(16)))/(T(18)*T(18))));
g1(16,18)=1-params(1)*params(13);
g1(16,19)=(-(T(19)*params(12)/(1+params(12))*1/params(12)));
g1(16,20)=(-(T(19)*1/params(36)));
g1(16,21)=(-T(19));
g1(17,2)=y(6)*y(8);
g1(17,6)=y(2)*y(8);
g1(17,7)=y(11);
g1(17,8)=y(6)*y(2);
g1(17,10)=(-1);
g1(17,11)=y(7);
g1(17,12)=1;
g1(18,18)=(-((1-params(13))*exp(y(18))*getPowerDeriv(exp(y(18)),T(20),1)*T(28)));
g1(18,19)=(-(T(24)*((-y(19))*(params(13)*T(21)*1/(y(19)*y(19))*log(y(21))+(1-params(13))*T(22)*1/(y(19)*y(19))*log(exp(y(18))))/T(23)-log(T(23)))));
g1(18,21)=1-T(28)*params(13)*getPowerDeriv(y(21),T(20),1);
g1(19,1)=1-params(17);
g1(20,19)=1/y(19)-params(18)*1/y(19);
g1(21,3)=1-params(6);
if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
end
end
