% 2016.10.22 GAL model v0.2

function dydt = GAL_sec(t,y,param)

% unwrap all parameters in param variable
% for fdn = fieldnames(param)'
%     eval( sprintf( '%s = param.%s;', fdn{1}, fdn{1} ) )
% end

% y(y<0) = 0;  % solve the non negative probelm

a1 = param.a1;
a2 = param.a2;
a3 = param.a3;
a4 = param.a4;
a80 = param.a80;
aR = param.aR;
ag1 = param.ag1;
ag2 = param.ag2;
ag3 = param.ag3;
ag4 = param.ag4;
ag80 = param.ag80;
d = param.d;

dsugar = param.dsugar;

KMglu = param.KMglu;
KMgal = param.KMgal;
% Kmt1 = param.Kmt1;
% Kmt2 = param.Kmt2;
kf3 = param.kf3;
kr3 = param.kr3;
kf83 = param.kf83;
kr83 = param.kr83;
kf84 = param.kf84;
kr84 = param.kr84;
kfR = param.kfR;
krR = param.krR;
kglu = param.kglu;
kgal = param.kgal;
% kct1 = param.kct1;
% kct2 = param.kct2;
KG1 = param.KG1;
KG2 = param.KG2;
KG3 = param.KG3;
KG80 = param.KG80;
KR1 = param.KR1;
KR3 = param.KR3;
KR4 = param.KR4;
n1 = param.n1;
n2 = param.n2;
n3 = param.n3;
n80 = param.n80;
nR1 = param.nR1;
nR3 = param.nR3;
nR4 = param.nR4;

exglu = param.exglu;
exgal = param.exgal;


% run ODE using un-wraped parameter
% y(1)=Gal1, y(2)=Gal2, y(3)=Gal3, y(4)=Gal4, y(5)=Gal80,
% y(6)=Gal3*,y(7)=Repressor,y(8)=R*,y(9)=C83,y(10)=C84,y(11)=glu,y(12)=gal;

dydt = zeros(12,1);
dydt(1) = a1+ag1*(y(4)^n1/(KG1^n1+y(4)^n1))*(KR1^nR1/(KR1^nR1+y(8)^nR1))-d*y(1);
dydt(2) = a2+ag2*(y(4)^n2/(KG2^n2+y(4)^n2))-d*y(2);
dydt(3) = a3+ag3*(y(4)^n3/(KG3^n3+y(4)^n3))*(KR3^nR3/(KR3^nR3+y(8)^nR3))-kf3*y(12)*y(3)+kr3*y(6)-d*y(3);
dydt(4) = a4+ag4*(KR4^nR4/(KR4^nR4+y(8)^nR4))-kf84*y(4)*y(5)+kr84*y(10)-d*y(4);
dydt(5) = a80+ag80*(y(4)^n80/(KG80^n80+y(4)^n80))-kf83*y(6)*y(5)+kr83*y(9)-kf84*y(4)*y(5)+kr84*y(10)-d*y(5);
dydt(6) = kf3*y(12)*y(3)-kr3*y(6)-kf83*y(6)*y(5)+kr83*y(9)-d*y(6);
dydt(7) = aR-kfR*y(11)*y(7)+krR*y(8)-d*y(7);
dydt(8) = kfR*y(11)*y(7)-krR*y(8)-d*y(8);
dydt(9) = kf83*y(6)*y(5)-kr83*y(9)-d*y(9);
dydt(10) = kf84*y(4)*y(5)-kr84*y(10)-d*y(10);
dydt(11) = kglu*(exglu/(KMglu/KMgal*exgal+exglu+KMglu))-dsugar*y(11);
dydt(12) = kgal*(exgal/(KMgal/KMglu*exglu+exgal+KMgal))-dsugar*y(12);

% if any(~isreal(dydt))
%   keyboard
% end