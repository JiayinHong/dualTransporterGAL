%% find out the approximate internal glucose level
% and appropriate value for KMglu

% the glucose gradient in one column fitting
% exglu = [56000000;28000000;14000000;7000000;3500000;1750000;875000;0];
exglu = linspace(0,10^8,20);    % test a broader range of glucose titration
kglu = 4350;    % parameter value from literature
dsugar = 7;

KM = 5.6*10^7;

figure
drawnow
for KMglu = KM * 0.5 .^ [0:1:5]     % try out a series of KMglu values
    inglu = kglu .* exglu ./ (exglu + KMglu) / dsugar;      % Michaelis Menten formula
    plot(exglu, inglu, 'LineWidth', 1.5)
    hold all
    pause(0.5)
end
set(gca, 'XTick', fliplr(exglu'))   % so that the value is increasing

% KMglu = 7*10^6;

%% decide appropriate value for KRs
nRs = 2;

figure
drawnow
for KRs = [200:100:1000]
    % the repression state mig1 fraction of total mig1
    Mig1s_frac = inglu .^ nRs ./ (KRs^nRs + inglu .^ nRs);

    plot(exglu, Mig1s_frac, 'LineWidth', 1.5)
    hold all
    pause(0.5)
end
set(gca, 'XTick', fliplr(exglu'))

% KRs = 300

%% glucose and galactose flux through GAL2
% the experiment setup for glucose gradient
exglu = [56000000;28000000;14000000;7000000;3500000;1750000;875000;0];
exgal = 28000000;

G2 = 1500;  % approximate level of GAL2
kG2 = 1;    % the scaling factor of GAL2 level

rG2 = 1;          % see 'The molecular genetics of hexose transport in yeasts'
KGglu = 2*10^6;   % Maier A, et al. (2002) 
KGgal = rG2*KGglu;  % Characteristics of Galactose Transport in Saccharomyces cerevisiae
                    % Cells and Reconstituted Lipid Vesicles
                    
glu_G2 = (kG2 * G2 .* exglu ./(1/rG2 * exgal + exglu + KGglu)) /7;  % dsugar = 7
gal_G2 = (kG2 * G2 * exgal ./(rG2 * exglu + exgal + KGgal)) /7;     % steady state sugar uptake through G2

%% glucose and galactose flux through HXT
kHXT = 4350;    % the maximum sugar uptake rate of HXT

rHXT = 200;     
KHXTglu = 7.5*10^6;
KHXTgal = rHXT*KHXTglu;

glu_hxt = (kHXT .* exglu ./(1/rHXT * exgal + exglu + KHXTglu)) /7;  % dsugar = 7
gal_hxt = (kHXT * exgal ./(rHXT * exglu + exgal + KHXTgal)) /7;     % steday state sugar uptake through HXT

% the complete expression of steady state intracellular sugar level
dglu = (kG2 * G2 .* exglu ./(1/rG2 * exgal + exglu + KGglu) + kHXT .* exglu ./(1/rHXT * exgal + exglu + KHXTglu)) /7;
dgal = (kG2 * G2 * exgal ./(rG2 .* exglu + exgal + KGgal) + kHXT * exgal ./(rHXT * exglu + exgal + KHXTgal)) /7;

%% HXT synthesis inhibited by G4
nHXT = 1;
G4 = logspace(0,1,20);
figure
drawnow
for KHXT = 3:3:21
    HXT_frac = KHXT .^nHXT ./ (KHXT .^nHXT + G4 .^nHXT);
    plot(G4,HXT_frac,'LineWidth', 1.5)
    hold all
    pause(0.3)
end

% KHXT = 9
