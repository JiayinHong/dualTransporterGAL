function stringout = changeunderscore( string, repstrinput )
%   This function is used to display underscore in titles normally,
%   instead of converting the letter next to the underscore a subscript

if iscell(string)
    stringChanged = {};
    for i = 1:length(string)
        stringChanged{i} = changeunderscore(string{i});
    end
    stringout = stringChanged;
    return
end

if nargin == 1;
    repstr = '\\\_';
else
    repstr = repstrinput;
end

stringout = regexprep(string, '_', repstr);

end

