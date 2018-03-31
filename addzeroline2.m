function [ hPlot ] = addzeroline2( varargin)
%addzeroline('flagdiag', 1, 'plotoption', {'linestyle', '-', 'color', 'k', 'linewidth', 1})
%addzeroline('plotoption', {'linestyle', '-', 'color', 'k', 'linewidth', 1})
% addzeroline('xpos', [], 'plotoption', {'linestyle', '-', 'linewidth', 1, 'color', 'k'})
%addzeroline Draw reference horizontal and vertical lines
%   Detailed explanation goes here
%   Written by Bo HUA,
%   Updated by Bo HUA, compatible update, add optional parameters xpos,
%   ypos, to define vertical and horizontal lines separately. Also re-write
%   the plotting part (variable definition), so that the new script will
%   be able to run compatibly with the old one.
%   Updated by Bo HUA, add default feature to always bury reference lines
%   to the back. This feature will be controlled by a flag variable.
%   Updated by Bo HUA, add hideline option, and left plot lines to the
%   back.

parser = inputParser;
addParamValue(parser,'figHandle', gca, @isnumeric);
addParamValue(parser,'pos',[0 0],@isnumeric);
addParamValue(parser,'xpos',nan,@isnumeric);
addParamValue(parser,'ypos',nan,@isnumeric);
addParamValue(parser,'hidelines',1,@isnumeric);
addParamValue(parser,'plotoption',{'color', 'k', 'linewidth', 1},@iscell);
addParamValue(parser,'flagDiag',0,@isnumeric);
addParamValue(parser,'flagDiagoffset',[],@isnumeric);

addParamValue(parser,'linecolor',rgb('gray'),@(x) isnumeric(x)|ischar(x));
addParamValue(parser,'linewidth',3,@isnumeric);
addParamValue(parser,'marker','none',@char);
addParamValue(parser,'linestyle','--',@char);

parse(parser, varargin{:});

pos = parser.Results.pos;
xpos = parser.Results.xpos;
ypos = parser.Results.ypos;
hidelines = parser.Results.hidelines;
plotoption = parser.Results.plotoption;
flagDiag = parser.Results.flagDiag;
flagDiagoffset = parser.Results.flagDiagoffset;

linecolor = parser.Results.linecolor;
linewidth = parser.Results.linewidth;
marker = parser.Results.marker;
linestyle = parser.Results.linestyle;

if flagDiag
    xpos = inf;
    ypos = inf;
end

if isnan(xpos) & isnan(ypos)
    xpos = pos(1);
    ypos = pos(2);
    
end
if ~isnan(xpos)
    if size(xpos, 2)>1 && size(xpos, 1)==1
        xpos = xpos';
    end
end

if ~isnan(ypos)
    
    if size(ypos, 2)>1 && size(ypos, 1)==1
        ypos = ypos';
    end
end

% plot(get(gca, 'xlim'), [pos(2) pos(2)], '--', 'linewidth', 3, 'color', rgb('gray'))
% plot([pos(1) pos(1)], get(gca, 'ylim'), '--', 'linewidth', 3, 'color', rgb('gray'))
hPlot = [];
xlimRange = get(gca, 'xlim');
ylimRange = get(gca, 'ylim');
% h1 = plot(xlimRange, repmat(pos(:,2),1,2)', '--', 'linewidth', 3, 'color', rgb('gray'));
% h2 = plot(repmat(pos(:,1),1,2)', ylimRange, '--', 'linewidth', 3, 'color', rgb('gray'));

if ~isempty(ypos) & ~isnan(ypos)
    h1 = plot(xlimRange, repmat(ypos,1,2)', 'linestyle', linestyle, 'marker', marker, 'linewidth', linewidth, 'color', linecolor,...
        plotoption{:});
    
    if hidelines
        
        tmp = get(gca, 'children');
        [~, a] = ismember(h1, tmp);
        tmp(a) = [];
        tmp = [tmp; h1];
        set(gca, 'children', tmp);
        
        %         tmp = tmp([(length(h1)+1):end, 1:length(h1)]);
        %         set(gca, 'children', tmp);
    end
else
    h1 = matlab.graphics.chart.primitive.Line.empty;
end

if ~isempty(xpos) & ~isnan(xpos)
    h2 = plot(repmat(xpos,1,2)', ylimRange, 'linestyle', linestyle, 'marker', marker,  'linewidth', linewidth, 'color', linecolor,...
        plotoption{:});
    if hidelines
        tmp = get(gca, 'children');
        [~, a] = ismember(h2, tmp);
        tmp(a) = [];
        tmp = [tmp; h2];
        set(gca, 'children', tmp);
    end
else
    h2 = matlab.graphics.chart.primitive.Line.empty;
end

set(gca, 'xlim', xlimRange)
set(gca, 'ylim', ylimRange)

hPlot = [h1;h2];

if flagDiag
    
    if isempty(flagDiagoffset)
        h3 = plot(get(gca, 'xlim'), get(gca, 'xlim'), 'linestyle', linestyle, 'marker', marker,  'linewidth', linewidth, 'color', linecolor,...
            plotoption{:});
        hPlot = [h1;h2;h3];
    else
        flagDiagoffset = flagDiagoffset(:);
        h3 = plot(repmat(get(gca, 'xlim'), length(flagDiagoffset), 1)',...
            repmat(get(gca, 'xlim'), length(flagDiagoffset), 1)' + repmat(flagDiagoffset, 1, 2)', 'linestyle', linestyle, 'marker', marker,  'linewidth', linewidth, 'color', linecolor,...
            plotoption{:});
        hPlot = [h1;h2;h3];
    end
end

end

