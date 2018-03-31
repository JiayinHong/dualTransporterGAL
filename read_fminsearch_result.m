function [param, obj, iter] = read_fminsearch_result(filepath)

fid = fopen(filepath);
if fid == -1
    error('file not exist')
end

% read the txt file, and use the new line to replace the old one
% until the end of the file
line = fgetl(fid);
line_minus_one = [];
while ~feof(fid)
    line_minus_two = line_minus_one;
    line_minus_one = line;
    line = fgetl(fid);
end

% convert the last line of the file into a one-row cell array,
% where each column is splited by tab
if strcmp(line, 'done')
    line = line_minus_two;
end

line = regexp( deblank(line), '\t', 'split' );
iter = str2double(line{1});
obj = str2double(line{2});
param = cellfun( @str2num, line(4:end) );

end
