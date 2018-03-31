function dy = test( t, y, a )
%   a simple exponential decay model
%   dy must be one column
%   or you can say dy = zeros(2,1)
%   then dy(1) = ...; dy(2) = ...;
dy(1,1) = -a * y(1);
dy(2,1) = y(1) * 100;

end
