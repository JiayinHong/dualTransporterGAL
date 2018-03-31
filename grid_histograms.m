yfprange = linspace(0,4,100);

galLabel = {'None','-8','-7','-6','-5','-4','-3','-2','-1','0','1','2'};
gluLabel = {'None','-6','-5','-4','-3','-2','-1','0'};

for iplate = 5
    [hf ha] = gridplot(8,12,100,100,'gaphorz',2,'gapvert',2);
    
    for r = 1:8
        for c = 1:12
            axes(ha(sub2ind([12 8],c,r)));
            dat = log10(alldata(iplate,9-r,c,2).yfp);
%             dat = log10(bc187_facs{9-r,c}.yfp);
            Y = histc(dat,yfprange)./numel(dat);
            plotfilledhist(yfprange,Y,'drawoutline',true);
            hold all

            xlim([0 4]);
            ylim([0 0.15]);
%             drawnow
            
            ax = gca;
            ax.XTickLabel = [];
            ax.YTickLabel = [];
            
            if r==8
                h1 = xlabel(galLabel{c});
                h1.FontSize = 12;
                h1.FontWeight = 'bold';
            end
            if c==1
                h2 = ylabel(gluLabel{9-r});
                h2.Rotation = 0;
                h2.HorizontalAlignment = 'right';
                h2.Margin = 10;
                h2.FontSize = 12;
                h2.FontWeight = 'bold';
            end
        end
    end
end

export_fig('WT_hist','-transparent','-c[NaN,NaN,NaN,NaN]')