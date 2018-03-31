function dydt = GALode3( t, y, param )
% 2017.08.03 GAL model v2.0
% updated from v1.0, only change previous parameter KMgal to alpha*KMglu,
% alpha is the ratio of KMgal over KMglu. This modification makes it easier
% to add a constrain on affinity constants of sugar binding to transporter.

% 2017.08.05 update
% use beta*kglu to replace kgal, beta is the ratio of kgal over kglu.
% Play the same trick on maximum transportation / uptake rate of sugar.

% 2017.08.08 GAL model v3.0
% Major update in regards of the formula of R and R*

y = max(y,0);  % solve the negative value problem
% this is more efficient than y(y<0)=0

% define variables
G1 = y(1);
G2 = y(2);
G3 = y(3);
G4 = y(4);
G80 = y(5);
G3s = y(6);
Mig1tot = y(7);
% Mig1s = y(8);
C83 = y(8);
C84 = y(9);
glu = y(10);
gal = y(11);

% unwrap parameters
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
alpha = param.alpha;
% KMgal = param.KMgal;
kf3 = param.kf3;
kr3 = param.kr3;
kf83 = param.kf83;
kr83 = param.kr83;
kf84 = param.kf84;
kr84 = param.kr84;
% kfR = param.kfR;
% krR = param.krR;
kglu = param.kglu;
beta = param.beta;
% kgal = param.kgal;
KG1 = param.KG1;
KG2 = param.KG2;
KG3 = param.KG3;
KG80 = param.KG80;
KR1 = param.KR1;
KR3 = param.KR3;
KR4 = param.KR4;
KRs = param.KRs;    % the affinity of glucose binding to free Mig1
n1 = param.n1;
n2 = param.n2;
n3 = param.n3;
n80 = param.n80;
nR1 = param.nR1;
nR3 = param.nR3;
nR4 = param.nR4;
nRs = param.nRs;    % the hill coefficient of glucose binding to free Mig1
exglu = param.exglu;
exgal = param.exgal;

Mig1s = glu^nRs/(KRs^nRs + glu^nRs) * Mig1tot;

G1_synthesis_term = G4^n1/(G4^n1 + KG1^n1) * KR1^nR1/(KR1^nR1 + Mig1s^nR1);
G2_synthesis_term = G4^n2/(G4^n2 + KG2^n2);
G3_synthesis_term = G4^n3/(G4^n3 + KG3^n3) * KR3^nR3/(KR3^nR3 + Mig1s^nR3);
G4_synthesis_term = KR4^nR4/(KR4^nR4 + Mig1s^nR4);
G80_synthesis_term = G4^n80/(G4^n80 + KG80^n80);

G3star_f = kf3 * gal * G3;
G3star_r = kr3 * G3s;

% Mig1star_f = kfR * glu * Mig1;
% Mig1star_r = krR * Mig1s;

C83_f = kf83 * G3s * G80;
C83_r = kr83 * C83;

C84_f = kf84 * G4 * G80;
C84_r = kr84 * C84;

dG1 = a1 + ag1 * G1_synthesis_term - d * G1;
dG2 = a2 + ag2 * G2_synthesis_term - d * G2;
dG3 = a3 + ag3 * G3_synthesis_term - G3star_f + G3star_r - d * G3;
dG4 = a4 + ag4 * G4_synthesis_term - C84_f + C84_r - d * G4;
dG80 = a80 + ag80 * G80_synthesis_term - C83_f + C83_r - C84_f + C84_r - d * G80;
dG3s = G3star_f - G3star_r - C83_f + C83_r - d * G3s;
dMig1tot = aR - d * Mig1tot;
% dMig1s = Mig1star_f - Mig1star_r - d * Mig1s;
dC83 = C83_f - C83_r - d * C83;
dC84 = C84_f - C84_r - d * C84;
% dglu = kglu * exglu/(1/alpha * exgal + exglu + KMglu) - dsugar * glu;
dglu = kglu * exglu/(exglu + KMglu) - dsugar * glu;
dgal = beta*kglu * exgal/(alpha * exglu + exgal + alpha*KMglu) - G3star_f + G3star_r  - dsugar * gal;

dydt = ...
    [ dG1 ...
    ; dG2 ...
    ; dG3 ...
    ; dG4 ...
    ; dG80 ...
    ; dG3s ...
    ; dMig1tot ...
    ; dC83 ...
    ; dC84 ...
    ; dglu ...
    ; dgal ];

end