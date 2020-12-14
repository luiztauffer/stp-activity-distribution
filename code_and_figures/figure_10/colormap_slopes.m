clear all; %close all;

N = 1; 
Nb = 1; 
Tb = 0.04;

for syntype = 1:3
    if syntype == 1
        %Facilitation --------
        tauf = .2;
        taud = .05;
        U = .1; 
    elseif syntype == 3
        %Depression --------
        tauf = .05;
        taud = .2;
        U = .7;
    elseif syntype == 2
        %Depression/Facilitation ----------
        tauf = 0.1;
        taud = 0.1;
        U = 0.4;   
    end
    
    %all_Rn = linspace(.5,10,20); 
    all_Rn = linspace(.5,10,20); 
    n = 0;
    for Rn = all_Rn
        n = n + 1;
        %all_max = Rn*linspace(.1,1,73);
        all_max = Rn*linspace(.1,1,73);
        k = 0;
        for re_max = all_max
            k = k + 1;
            all_Re = linspace(0,re_max,10);
            nCols = length(all_Re);
            NRV = zeros(1,nCols);

            j = 1;
            NRVbasal = solve_NRV(tauf, taud, U, Rn, 0, N, 1, Tb);
            for Re = all_Re
                 NRV(j) = solve_NRV(tauf, taud, U, Rn, Re, N, Nb, Tb) - NRVbasal;
                j = j+1;
            end


            [p,~,~,~,stats] = regress(NRV',[ones(length(all_Re),1) all_Re']);
            slopes(n,k) = p(2);
            rsqu(n,k) = 100*stats(1);
            
            
            if k == 1
                slope_ref = p(2);
            end

        end
        
        %Absolut value of the Slope Deviation
        sl_dev(n,:) = abs( 100*(slopes(n,:)-slope_ref)/slope_ref );
        
        %Find where deviation = 1%
        [~,ind1perc(n,syntype)] = min( abs(sl_dev(n,:)-1) );
    end
    
    nCols = length(all_max);
    nRows = length(all_Rn);

    
    if syntype == 1
        sb(4) = subplot('Position',[.08 .4 .22 .23]);
    elseif syntype == 2
        sb(5) = subplot('Position',[.39 .4 .22 .23]);
    elseif syntype == 3
        sb(6) = subplot('Position',[.7 .4 .22 .23]);
    end
    imagesc(sl_dev);
    set(gca,'YDir','normal'); hold on;
    plot(ind1perc(:,syntype),1:nRows,'color',[.7,.7,.7],'Linewidth',2);
    colormap(sb(3+syntype),'jet'); 
    h = gca();
    h.XTick = linspace(1,nCols,4);
    h.YTick = linspace(1,nRows,5);
    h.XTickLabel = strread(num2str( [.1,.4,.7,1] ),'%s');
    h.YTickLabel = strread(num2str( all_Rn(round(linspace(1,nRows,5))) ),'%s');
    xlabel('r_{\delta}/r_{bas}');
    ylabel('r_{bas} [Hz]');
    
    
    if syntype == 1
        sb(7) = subplot('Position',[.08 .08 .22 .23]);
    elseif syntype == 2
        sb(8) = subplot('Position',[.39 .08 .22 .23]);
    elseif syntype == 3
        sb(9) = subplot('Position',[.7 .08 .22 .23]);
    end
    imagesc(rsqu);
    set(gca,'YDir','normal'); hold on;
    colormap(sb(6+syntype),'bone'); 
    h = gca();
    h.XTick = linspace(1,nCols,4);
    h.YTick = linspace(1,nRows,5);
    h.XTickLabel = strread(num2str( [.1,.4,.7,1] ),'%s');
    h.YTickLabel = strread(num2str( all_Rn(round(linspace(1,nRows,5))) ),'%s');
    xlabel('r_{\delta}/r_{bas}');
    ylabel('r_{bas} [Hz]');
    
    
    allmin1(syntype) = min(min(sl_dev));
    allmax1(syntype) = max(max(sl_dev));
    allmin2(syntype) = min(min(rsqu));
    allmax2(syntype) = max(max(rsqu));
end

minc1 = min(allmin1);
maxc1 = max(allmax1);
minc2 = min(allmin2);
maxc2 = max(allmax2);
for si = 4:6
    set(sb(si),'clim',[minc1,maxc1]);
end
for si = 7:9
    set(sb(si),'clim',[minc2,maxc2]);
end

c1 = colorbar(sb(4));
set(c1,'Position',[.925 .4 .02 .23]);
title(c1,'|S^s dev| [%]');

c2 = colorbar(sb(7));
set(c2,'Position',[.925 .08 .02 .23]);
set(c2,'YTick',[99.9, maxc2]);
set(c2,'XTickLabel',{'99.9', '100'});
title(c2,'R-squared [%]');

