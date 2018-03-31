function stop = fmincon_outfun(x, optimValues, state, filepath)
% user-defined function that the optimization function calls at each
% iteration, for more details, see help documentation for 'Output
% Function', 'Optimization Options Reference', 'optimset' and 'optimoptions'
stop = false;

switch state
    case 'init'
        fid = fopen(filepath, 'w');
    otherwise
        fid = fopen(filepath, 'a+');
end

fprintf(fid, '%d\t', optimValues.iteration);
fprintf(fid, '%1.3e\t', optimValues.fval);
fprintf(fid, '%d\t', optimValues.funccount);
fprintf(fid, '%1.3e\t', optimValues.gradient);
fprintf(fid, '%1.3e\t', optimValues.stepsize);
% funccount is cumulative number of function evaluations
% fval is function value at current point.
% gradient is current gradient of objective function

for i = 1:length(x)
    fprintf(fid, '%1.3e\t', x(i));
end

fprintf(fid, '\n');

switch state
    case 'done'
        fprintf(fid, 'done\n');
end

fclose(fid);