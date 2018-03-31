function stop = fminsearch_outfun(x, optimValues, state, filepath)

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

for i = 1:length(x)
    fprintf(fid, '%1.3e\t', x(i));
end

fprintf(fid, '\n');

switch state
    case 'done'
        fprintf(fid, '\ndone\n');
end

fclose(fid);