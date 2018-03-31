function param = set_parameter(param_set_i)

param = struct();

switch param_set_i
    case 1      % parameter set I is inferred from reference
        % the following production rate and dilution rate are estimated
        % from experimental data
        %         param.a1 = 0;
        param.a1 = 0.001; % not able to generate a random number from lognormal distribution whose mu=log0
        param.a2 = 0;
        param.a3 = 0.25 * 2;
        param.a4 = 0.058;
        param.a80 = 0.272;
        param.aR = 0.288;
        param.ag1 = 51.77;
        param.ag2 = 24.36;      % no experimental data available
        param.ag3 = 2.045;
        param.ag4 = 0.011;
        param.ag80 = 0.735;
        param.d = 0.0077;
        %         param.dsugar = 0.077;   % no experimental data available
        param.dsugar = 7;
        
        % the following association rate and dissociation rate referred to
        % Venturelli's model
        param.kf3 = 70.3;
        param.kr3 = 3391;
        param.kf83 = 41.1 * 1000;
        param.kr83 = 700.1;
        param.kf84 = 95.2 * 10;
        param.kr84 = 1237;
        param.kfR = 56.2;
        param.krR = 3564;
        
        % the following transcription KMs referred to Venturelli's model
        param.KG1 = 41.6;
        param.KG2 = 66.6;       % no experimental data available
        param.KG3 = 32.2;       % no experimental data available
        param.KG80 = 14;
        param.KR1 = 67.4;
        param.KR3 = 18.8;       % no experimental data available
        param.KR4 = 33.8;
        
        % the following sugar transportation rate and KMs referred to
        % Bennett's model, the KMs were adapted based on experimental data
        param.kglu = 4350;
        param.kgal = 4350;
        param.KMglu = 7.5*10^5;
        param.KMgal = 3*10^6 * 10;  % Aug1st update, more biological plausible
        % param.KM1 = 7.5*10^5;
        % param.KM2 = 3*10^6;
        
        % the following Hill coefficient referred to Venturelli's model
        param.n1 = 3;
        param.n2 = 3;       % no experimental data available
        param.n3 = 2;
        param.n80 = 2;
        param.nR1 = 2;
        param.nR3 = 2;      % no experimental data available
        param.nR4 = 1;
        
    case 2      % the current default values, manually tuned from set I.
        param.a1 = 0.001;
        param.a2 = 0;
        param.a3 = 0.15;
        param.a4 = 0.058;
        param.a80 = 0.272;
        param.aR = 0.288;
        param.ag1 = 10;
        param.ag2 = 24.36;
        param.ag3 = 1.045;
        param.ag4 = 0.011;
        param.ag80 = 0.735;
        param.d = 0.0077;
        param.dsugar = 70;
        param.kf3 = 2.4;
        param.kr3 = 40000;
        param.kf83 = 720;
        param.kr83 = 60;
        param.kf84 = 285.6;
        param.kr84 = 1237;
        % param.kfR = 0.0562;
        % param.krR = 1069.2;
        param.KG1 = 3.67;
        param.KG2 = 66.6;
        param.KG3 = 6.67;
        param.KG80 = 24.03;
        param.KR1 = 17.58;
        param.KR3 = 28.8;
        param.KR4 = 6.74;
        param.KRs = 100;
        param.kglu = 4350;
        param.beta = 1;
        % param.kgal = 4350;
        param.KMglu = 7.5*10^5;
        param.alpha = 200;
        % param.KMgal = 3*10^6;
        param.n1 = 3;
        param.n2 = 3;
        param.n3 = 2;
        param.n80 = 2;
        param.nR1 = 2;
        param.nR3 = 2;
        param.nR4 = 1;
        param.nRs = 2;
        
    case 3      % use alpha*KMglu to replace KMgal
        % use beta*kglu to replace kgal
        % remove kfR, krR; add KRs, nRs
        % put back R*, internal gluc level not depend on competitive uptake
        % but on external gluc alone
        
        param.a1 = 0.001; % not able to generate a random number from lognormal distribution whose mu=log0
        param.a2 = 0;
        param.a3 = 0.25 * 2;
        param.a4 = 0.058;
        param.a80 = 0.272;
        param.aR = 0.288;
        param.ag1 = 51.77;
        param.ag2 = 24.36;      % no experimental data available
        param.ag3 = 2.045;
        param.ag4 = 0.011;
        param.ag80 = 0.735;
        param.d = 0.0077;
        param.dsugar = 7;
        
        % the following association rate and dissociation rate referred to
        % Venturelli's model
        param.kf3 = 70.3;
        param.kr3 = 3391;
        param.kf83 = 41.1 * 1000;
        param.kr83 = 700.1;
        param.kf84 = 95.2 * 10;
        param.kr84 = 1237;
        
        % the following transcription KMs referred to Venturelli's model
        param.KG1 = 41.6;
        param.KG2 = 66.6;       % no experimental data available
        param.KG3 = 32.2;       % no experimental data available
        param.KG80 = 14;
        param.KR1 = 67.4;
        param.KR3 = 18.8;       % no experimental data available
        param.KR4 = 33.8;
        param.KRs = 300;        % comparable value to internal gluc
        
        % the following sugar transportation rate and KMs referred to
        % Bennett's model, the KMs were adapted based on experimental data
        param.kglu = 4350;
        param.beta = 1;
        %         param.kgal = 4350;
        param.KMglu = 7.5*10^5 * 10;    % to get rid of eternal saturation
        param.alpha = 200;
        
        % the following Hill coefficient referred to Venturelli's model
        param.n1 = 3;
        param.n2 = 3;       % no experimental data available
        param.n3 = 2;
        param.n80 = 2;
        param.nR1 = 1;
        param.nR3 = 1;      % no experimental data available
        param.nR4 = 1;
        param.nRs = 2;
        
    case 4      % use alpha*KMglu to replace KMgal
        % use beta*kglu to replace kgal
        
        param.a1 = 0.001; % not able to generate a random number from lognormal distribution whose mu=log0
        param.a2 = 0;
        param.a3 = 0.25 * 2;
        param.a4 = 0.058;
        param.a80 = 0.272;
        param.aR = 0.288;
        param.ag1 = 51.77;
        param.ag2 = 24.36;      % no experimental data available
        param.ag3 = 2.045;
        param.ag4 = 0.011;
        param.ag80 = 0.735;
        param.d = 0.0077;
        param.dsugar = 7;
        
        % the following association rate and dissociation rate referred to
        % Venturelli's model
        param.kf3 = 70.3;
        param.kr3 = 3391;
        param.kf83 = 41.1 * 1000;
        param.kr83 = 700.1;
        param.kf84 = 95.2 * 10;
        param.kr84 = 1237;
        param.kfR = 56.2;
        param.krR = 3564;
        
        % the following transcription KMs referred to Venturelli's model
        param.KG1 = 41.6;
        param.KG2 = 66.6;       % no experimental data available
        param.KG3 = 32.2;       % no experimental data available
        param.KG80 = 14;
        param.KR1 = 67.4;
        param.KR3 = 18.8;       % no experimental data available
        param.KR4 = 33.8;
        
        % the following sugar transportation rate and KMs referred to
        % Bennett's model, the KMs were adapted based on experimental data
        param.kglu = 4350;
        param.beta = 1;
        %         param.kgal = 4350;
        param.KMglu = 7.5*10^5;
        param.alpha = 200;
        %         param.KMgal = 3*10^6 * 10;  % Aug1st update, more biological plausible
        
        % the following Hill coefficient referred to Venturelli's model
        param.n1 = 3;
        param.n2 = 3;       % no experimental data available
        param.n3 = 2;
        param.n80 = 2;
        param.nR1 = 2;
        param.nR3 = 2;      % no experimental data available
        param.nR4 = 1;
        
    case 5      % use alpha*KMglu to replace KMgal
        % use beta*kglu to replace kgal
        % remove kfR, krR; add KRs, nRs
        
        param.a1 = 0.001; % not able to generate a random number from lognormal distribution whose mu=log0
        param.a2 = 0;
        param.a3 = 0.25 * 2;
        param.a4 = 0.058;
        param.a80 = 0.272;
        param.aR = 0.288;
        param.ag1 = 51.77;
        param.ag2 = 24.36;      % no experimental data available
        param.ag3 = 2.045;
        param.ag4 = 0.011;
        param.ag80 = 0.735;
        param.d = 0.0077;
        param.dsugar = 7;
        
        % the following association rate and dissociation rate referred to
        % Venturelli's model
        param.kf3 = 70.3;
        param.kr3 = 3391;
        param.kf83 = 41.1 * 1000;
        param.kr83 = 700.1;
        param.kf84 = 95.2 * 10;
        param.kr84 = 1237;
        %         param.kfR = 56.2;
        %         param.krR = 3564;
        
        % the following transcription KMs referred to Venturelli's model
        param.KG1 = 41.6;
        param.KG2 = 66.6;       % no experimental data available
        param.KG3 = 32.2;       % no experimental data available
        param.KG80 = 14;
        param.KR1 = 67.4;
        param.KR3 = 18.8;       % no experimental data available
        param.KR4 = 33.8;
        param.KRs = 10;        % arbitrary picked value
        
        % the following sugar transportation rate and KMs referred to
        % Bennett's model, the KMs were adapted based on experimental data
        param.kglu = 4350;
        param.beta = 1;
        %         param.kgal = 4350;
        param.KMglu = 7.5*10^5;
        param.alpha = 200;
        %         param.KMgal = 3*10^6 * 10;  % Aug1st update, more biological plausible
        
        % the following Hill coefficient referred to Venturelli's model
        param.n1 = 3;
        param.n2 = 3;       % no experimental data available
        param.n3 = 2;
        param.n80 = 2;
        param.nR1 = 2;
        param.nR3 = 2;      % no experimental data available
        param.nR4 = 1;
        param.nRs = 2;
        
    case 6      % use alpha*KMglu to replace KMgal
        % use beta*kglu to replace kgal
        % remove kfR, krR; add KRs, nRs
        % make nR1, nR3, and nR4 = 1
        
        param.a1 = 0.001; % not able to generate a random number from lognormal distribution whose mu=log0
        param.a2 = 0;
        param.a3 = 0.25 * 2;
        param.a4 = 0.058;
        param.a80 = 0.272;
        param.aR = 0.288;
        param.ag1 = 51.77;
        param.ag2 = 24.36;      % no experimental data available
        param.ag3 = 2.045;
        param.ag4 = 0.011;
        param.ag80 = 0.735;
        param.d = 0.0077;
        param.dsugar = 7;
        
        % the following association rate and dissociation rate referred to
        % Venturelli's model
        param.kf3 = 70.3;
        param.kr3 = 3391;
        param.kf83 = 41.1 * 1000;
        param.kr83 = 700.1;
        param.kf84 = 95.2 * 10;
        param.kr84 = 1237;
        
        % the following transcription KMs referred to Venturelli's model
        param.KG1 = 41.6;
        param.KG2 = 66.6;       % no experimental data available
        param.KG3 = 32.2;       % no experimental data available
        param.KG80 = 14;
        param.KR1 = 67.4 * 10^4;
        param.KR3 = 18.8 * 10^4;       % no experimental data available
        param.KR4 = 33.8 * 10^4;
        param.KRs = 10;        % arbitrary picked value
        
        % the following sugar transportation rate and KMs referred to
        % Bennett's model, the KMs were adapted based on experimental data
        param.kglu = 4350;
        param.beta = 1;
        param.KMglu = 7.5*10^5;
        param.alpha = 200;
        %         param.KMgal = 3*10^6 * 10;  % Aug1st update, more biological plausible
        
        % the following Hill coefficient referred to Venturelli's model
        param.n1 = 3;
        param.n2 = 3;       % no experimental data available
        param.n3 = 2;
        param.n80 = 2;
        param.nR1 = 1;
        param.nR3 = 1;      % no experimental data available
        param.nR4 = 1;
        param.nRs = 2;
        
    case 7      % major change - two transporters version
        
        param.a1 = 0.001; % not able to generate a random number from lognormal distribution whose mu=log0
        param.a2 = 0;
        param.a3 = 0.25 * 2;
        param.a4 = 0.058;
        param.a80 = 0.272;
        param.aR = 0.288;
        param.ag1 = 51.77;
        param.ag2 = 24.36;      % no experimental data available
        param.ag3 = 2.045;
        param.ag4 = 0.011;
        param.ag80 = 0.735;
        param.d = 0.0077;
        param.dsugar = 7;
        
        % the following association rate and dissociation rate referred to
        % Venturelli's model
        param.kf3 = 70.3;
        param.kr3 = 3391;
        param.kf83 = 41.1 * 1000;
        param.kr83 = 700.1;
        param.kf84 = 95.2 * 10;
        param.kr84 = 1237;
        
        % the following transcription KMs referred to Venturelli's model
        param.KG1 = 41.6;
        param.KG2 = 66.6;       % no experimental data available
        param.KG3 = 32.2;       % no experimental data available
        param.KG80 = 14;
        param.KR1 = 67.4;
        param.KR3 = 18.8;       % no experimental data available
        param.KR4 = 33.8;
        param.KRs = 300;        % comparable value to internal gluc
        
        % the following new parameters are deduced from various literature
        param.kG2 = 1;      % the scaling factor of GAL2 level
        param.kHXT = 4350;  % steal the value from previous 'kglu'
        param.rG2 = 1;      % presume GAL2 has the same affinity on binding glucose and galactose
        param.rHXT = 200;   % steal the value from previous 'alpha'
        param.KGglu = 2*10^6;   % see 'The molecular genetics of hexose transport in yeasts'
        param.KHXTglu = 7.5*10^6;   % steal the value from previous 'KMglu'
        
        % the following Hill coefficient referred to Venturelli's model
        param.n1 = 3;
        param.n2 = 2;       % no experimental data available
        param.n3 = 2;
        param.n80 = 2;
        param.nR1 = 1;
        param.nR3 = 1;      % no experimental data available
        param.nR4 = 1;
        param.nRs = 2;
        
    case 8      % major change - two transporters version + HXT species
                % compatible with GALode5
        
        param.a1 = 0.001; % not able to generate a random number from lognormal distribution whose mu=log0
        param.a2 = 0.08;
        param.a3 = 0.25;    % used to be 0.5
        param.a4 = 0.063;
        param.a80 = 0.282;
        param.aR = 0.295;
        param.ag1 = 350;    % param.ag1 = 51.77;
        param.ag2 = 11.28;  % param.ag2 = 24.36;  no experimental data available
        param.ag3 = 1.92;   % param.ag3 = 2.045;
        param.ag4 = 0.007;
        param.ag80 = 0.389; % param.ag80 = 0.735;
        
        param.a0HXT = 0.06;
        param.aHXT = 11.32;
        
        param.d = 0.0089;   % 0.0077 inferred from literature
        param.dsugar = 7;
        
        % the following association rate and dissociation rate referred to
        % Venturelli's model
        param.kf3 = 75.16;
        param.kr3 = 3391;
        param.kf83 = 52656; % param.kf83 = 41100;
        param.kr83 = 700.1;
        param.kf84 = 502;   % param.kf84 = 952;
        param.kr84 = 1237;
        
        % the following transcription KMs referred to Venturelli's model
        param.KG1 = 16.82;
        param.KG2 = 23.51;  % param.KG2 = 66.6;  no experimental data available
        param.KG3 = 17.3;   % param.KG3 = 32.2;  no experimental data available
        param.KG80 = 12.92;
        param.KR1 = 2.83;   % param.KR1 = 67.4;
        param.KR3 = 5.45;       % no experimental data available
        param.KR4 = 15.76;
        param.KRs = 31.34;  % param.KRs = 300;  comparable value to internal gluc
        
        param.KHXT = 5.09;
        
        % the following new parameters are deduced from various literature
        param.kG2 = 1.2;      % the scaling factor of GAL2 level
        
        param.rcat = 1;     % presum GAL2 and HXT has the same maximum uptake rate
        param.rG2 = 1;      % presume GAL2 has the same affinity on binding glucose and galactose
        param.rHXT = 177;   % param.rHXT = 200; steal the value from previous 'alpha'
        param.KGglu = 2*10^6;   % see 'The molecular genetics of hexose transport in yeasts'
        param.KHXTglu = 2*10^6; % param.KHXTglu = 7.5*10^6; steal the value from previous 'KMglu'
        
        % the following Hill coefficient referred to Venturelli's model
        param.n1 = 3;
        param.n2 = 2;       % no experimental data available
        param.n3 = 2;
        param.n80 = 2;
        param.nR1 = 2;
        param.nR3 = 1;      % no experimental data available
        param.nR4 = 1;
        param.nRs = 2;
        
        param.nHXT = 10;    % though strange, based on MCMC search. (prior set to be 2)
end
