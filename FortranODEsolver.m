function [ y ] = FortranODEsolver( param, y0, t )
%   This function first unwrap all the parameters, then input them to
%   Fortran ODE solver.
%   2017.07.31 by JH

p = MiniDataFrame(param);   % convert parameter struct to MiniDataFrame

p.KMgal = p.KMglu .* p.alpha;
p.kgal = p.kglu .* p.beta;
% p = [a1, a2, a3, a4, a80, aR, ag1, ag2, ag3, ag4, ag80 ...
%     , d, dsugar, KMglu, KMgal, kf3, kr3, kf83, kr83, kf84, kr84 ...
%     , kfR, krR, kglu, kgal, KG1, KG2, KG3, KG80, KR1, KR3, KR4 ...
%     , n1, n2, n3, n80, nR1, nR3, nR4, exglu, exgal];

% call the function 'galsim' to solve the ode
y = galsim(p, y0, t);
    
% get the steady state value of each species
% yt = y(end, :);

end
