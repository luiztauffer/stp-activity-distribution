%This script plots the Figure 2
clear all; %close all;

%Part 1 - Single fiber ----------------------------------------------------
figure(); set(gcf,'color','w','Position', [50, 50, 1000, 700]);
for syntype = 1:2
    if syntype == 2
        %Depression --------
        tauf = .05;
        taud = .2;
        U = .7;
    elseif syntype == 1
        %Facilitation --------
        tauf = .2;
        taud = .05;
        U = .1;    
    end

    Rn = 0.5; 
    N = 1; 
    Nb = 1; 
    Tb = 0.04;

    all_Re = linspace(50,200,3);
    
    dt = 0.0001;
    xx = 1000*(dt:dt:2*Tb);
    sb0(syntype) = subplot(4,4,4*syntype-3);
    ic = 1;
    for Re = all_Re
        [y, upb, Rseq] = solve_NRV_seq(tauf, taud, U, Rn, Re, N, 1, Tb);
        PRR = y(2,:).*upb.*Rseq;
        PRR(1) = [];
        plot(xx,PRR,'color',1-ic*[.25,.25,.25],'Linewidth',2); hold on;
        ic = ic+1;
    end
    plot([20,20],[0,1.1*max(PRR)],'--k');
    plot([60,60],[0,1.1*max(PRR)],'--k');
    ylabel('PRR^{s}(t)/sec'); box off;
    ylim([0,1.1*max(PRR)]); xlim([0,80]);
    h = gca(); h.XTick = [20, 60];
    if syntype == 1
        l1 = legend('50 Hz','125 Hz','200 Hz'); legend boxoff;
        set(l1,'Position',[0.16 0.89 0.1 0.016]);
        title('Facilitation dominated','FontWeight','Normal');
        set(gca,'Position',[.07 .81 .12 .13]);
    else
        xlabel('Time [ms]');  
        title('Depression dominated','FontWeight','Normal');
        set(gca,'Position',[.07 .6 .12 .13]);
    end
    
    
    all_Re = linspace(1,500,500);       
    nCols = length(all_Re);
    NRV = zeros(1,nCols);
    
    j = 1;
    NRVbasal = solve_NRV(tauf, taud, U, Rn, 0, N, 1, Tb);
    for Re = all_Re
        NRV(j) = solve_NRV(tauf, taud, U, Rn, Re, N, Nb, Tb) - NRVbasal;
        j = j+1;
    end
    
    sb3 = subplot('Position',[.32 .6 .18 .35]);
    plot(all_Re,NRV','Linewidth',2); hold on;
    if syntype == 2
        %plot(all_Re,all_Re*Tb,'k','Linewidth',2);
        xlim([0, 500]); ylim([0, 1.6]); box off;
        l2 = legend('Facilitation','Depression'); legend boxoff;
        set(l2,'Position',[0.38 0.67 0.13 0.059]);
        xlabel('r_{ext} [Hz]');
        ylabel('Q^{s}_{ext}');
    end
    
    [optfreq, optgain] = theoretical_optfreq(tauf, taud, U, Rn, Tb);
    subplot('Position',[.75 .1 .2 .14]);
    semilogx(linspace(0.01,1,1000),optgain,'Linewidth',2); hold on; 
    box off; %ylim([0, 140]); 
    h = gca(); h.XTickLabel = strread(num2str([1 10 100]),'%s');
    xlabel('r_{\delta} [% of r_{bas}]');
    ylabel('G_{opt} [%]');
    ylim([-2,100]);
end


%Part 2 - Acumulated population PRR ---------------------------------------
tauf = .2;
taud = .05;
U = .1;    

Rn = .5; 
N = 160000; 
Tb = 0.04;

Re = N*Rn*0.08;
all_Nb = unique( logspace(0,5,10000) );

nCols = length(all_Nb);
NRV = zeros(1,nCols);
NRV2 = zeros(1,nCols);

i = 0;
for Nb = all_Nb
    i = i+1;
    NRV(i) = solve_NRV(tauf, taud, U, Rn, Re, N, Nb, Tb);
    NRV2(i) = solve_NRV(tauf, taud, U, Rn, 2*Re, N, Nb, Tb);
end

NRVbasal = solve_NRV(tauf, taud, U, Rn, 0, N, 1, Tb);
NRV = NRV - NRVbasal;
NRV2 = NRV2 - NRVbasal;
[mval, mind] = max(NRV);


sb4 = subplot('Position',[.6 .6 .3 .36]);
semilogx(all_Nb,NRV,'k','LineWidth',2); hold on;
semilogx(all_Nb,NRV(end)*ones(1,length(all_Nb)),'--k'); 
semilogx(all_Nb,NRV2,'k','LineWidth',2); hold on;
semilogx(all_Nb,NRV2(end)*ones(1,length(all_Nb)),'--k'); box off;
text(300,40,'R_{ext}=8 kHz');
text(700,80,'R_{ext}=16 kHz');
%hmax = plot(all_Nb(mind), mval, 'ok','LineWidth',2);
%l00 = legend(hmax,['N_{opt} = ',num2str(round(all_Nb(mind)))]); legend boxoff;
%set(l00,'Position',[0.75 0.9 0.15 0.05]);
xticks([10, 1000, 100000]);
xlabel('N_{ext}'); ylabel('Q^p_{ext}-Q^p_{bas}');



set(sb0(1),'Position',[.08 .81 .16 .13]);
set(sb0(2),'Position',[.08 .6 .16 .13]);
set(sb3,'Position',[.39 .6 .21 .35]);
set(sb4,'Position',[.70 .6 .23 .36]);

set(l1,'FontSize',12,'Position',[0.2 0.89 0.1 0.016]);
set(l2,'FontSize',12,'Position',[0.46 0.67 0.13 0.059]);
%set(l00,'FontSize',12,'Position',[0.8 0.895 0.08 0.06]);



%Part 3 - Lower figures ---------------------------------------
mapping_NRV;


fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',12);
set(findall(fig,'-property','FontName'),'FontName','Arial');