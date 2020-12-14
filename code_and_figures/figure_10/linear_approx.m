%Single fiber Q(r) linear approximation for small changes in r
clear all; %close all;

Rn = .5; 
N = 1; 
Nb = 1; 
Tb = 0.04;

all_max = Rn*linspace(.1,1,30);%*logspace(-1,1,100);  

figure(); set(gcf,'color','w','Position', [50, 0, 1000, 700]);
k = 0;
for re_max = all_max
    k = k + 1;
    for syntype = 1:3
        if syntype == 1
            %Facilitation --------
            tauf = .2;
            taud = .05;
            U = .1; 
        elseif syntype == 2
            %Depression --------
            tauf = .05;
            taud = .2;
            U = .7;
        elseif syntype == 3
            %Depression/Facilitation ----------
            tauf = 0.1;
            taud = 0.1;
            U = 0.4;   
        end

        all_Re = linspace(0,re_max,101);
        nCols = length(all_Re);
        NRV = zeros(1,nCols);

        j = 1;
        NRVbasal = solve_NRV(tauf, taud, U, Rn, 0, N, 1, Tb);
        for Re = all_Re
             NRV(j) = solve_NRV(tauf, taud, U, Rn, Re, N, Nb, Tb) - NRVbasal;
            j = j+1;
        end

        %[p,S] = polyfit(all_Re,NRV,1);
        %slopes(k,syntype) = p(1);
        %[y,delta] = polyval(p,all_Re,S);
        %errors(k,syntype) = mean(delta);
        [p,~,~,~,stats] = regress(NRV',[ones(length(all_Re),1) all_Re']);
        slopes(k,syntype) = p(2);
        rsqu(k,syntype) = 100*stats(1);
        
        if k == 1
            slope_ref(syntype) = p(2);
        end
        
        if re_max == max(all_max)
            sb(1) = subplot('Position',[.08 .72 .22 .23]);
            plot(all_Re/Rn,NRV','Linewidth',2); hold on;
            %plot(all_Re,y,'Linewidth',.8,'color',[.5,.5,.5]);
            %plot(all_Re,y+delta,'Linewidth',.8,'color',[.5,.5,.5]);
            %plot(all_Re,y-delta,'Linewidth',.8,'color',[.5,.5,.5]);
            if syntype == 3
                box off;
                l1 = legend('Fac','Dep','Fac/Dep'); legend boxoff;
                set(l1,'Position',[0.085 0.85 0.13 0.1]);
                sb(1).XTick = [.1, .4, .7, 1];
                xlabel('r_{\delta}/r_{bas}');
                ylabel('Q^s_{\delta}-Q^s_{bas}');
            end
        end
        
    end
    
    
end
 
sb(2) = subplot('Position',[.39 .72 .22 .23]);
plot(all_max/Rn,100*(slopes(:,1)-slope_ref(1))/slope_ref(1),'Linewidth',2); hold on;
plot(all_max/Rn,100*(slopes(:,2)-slope_ref(2))/slope_ref(2),'Linewidth',2);
plot(all_max/Rn,100*(slopes(:,3)-slope_ref(3))/slope_ref(3),'Linewidth',2); box off;
sb(2).XTick = [.1, .4, .7, 1];
ylim([-.6,.6]);
ylabel('S^s deviation [%]');
xlabel('r_{\delta}/r_{bas}');


sb(3) = subplot('Position',[.72 .72 .22 .23]);
plot(all_max/Rn,rsqu,'Linewidth',2); box off;
ylim([99.998, 100]);
sb(3).XTick = [.1, .4, .7, 1];
sb(3).YTick = [99.998, 99.999, 100];
xlabel('r_{\delta}/r_{bas}');
ylabel('R-squared [%]');

colormap_slopes;

fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',12);
set(findall(fig,'-property','FontName'),'FontName','Arial');






