function stop = bads_outfun(x, optimState, state, filepath)

stop = false;

switch state
    case 'init'
        fid = fopen(filepath, 'w');
    otherwise
        fid = fopen(filepath, 'a+');
end

fprintf(fid, '%.1fh\t', optimState.totalfunevaltime/3600);
fprintf(fid, '%d\t', optimState.iter);
fprintf(fid, '%.3f\t', optimState.fval);
fprintf(fid, '%d\t', optimState.funccount);

for i = 1:length(x)
    fprintf(fid, '%1.3e\t', x(i));
end

fprintf(fid, '\n');

switch state
    case 'done'
        fprintf(fid, '\nThe bads optimization search is done!');
end

fclose(fid);