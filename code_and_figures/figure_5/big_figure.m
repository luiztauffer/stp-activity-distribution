clear all; close all;

figure(); set(gcf,'color','w','Position', [50, 35, 1000, 750]);

mymap = my_cmap(-90,120,0);   %Blue to Red, with center on 0

for syntype = 1:2
    if syntype == 1  %Input to target --------
        tauf = .2;
        taud = .05;
        U = 0.1; 
        syntitle = 's1 - facilitatory';        
    elseif syntype == 2  %Input to Inhibitory neuron --------
        tauf = .05;
        taud = .2;
        U = .7;
        syntitle = 's2 - depressing';
    end

    % Cerebellum - Parallel fibers -> Purkinje cells
    Rn = .5; 
    N = 160000; 
    Tb = 0.04;
    
%     % Hipoccampus - CA3 -> CA1
%     Rn = 5.; 
%     N = 1000; 
%     Tb = 0.200;


    all_Nb = unique( logspace(0,4,1000) );
    %all_Re = linspace(1000,8000,71);
    all_Re = linspace(1600,8000,160);

    
    NRVbasal = solve_NRV(tauf, taud, U, Rn, 0, N, 1, Tb);
    
    nRows = length(all_Nb);
    nCols = length(all_Re);
    NRV = zeros(nRows,nCols);
    i = 1;
    for Nb = all_Nb
        j = 1;
        for Re = all_Re
            NRV(i,j) = solve_NRV(tauf, taud, U, Rn, Re, N, Nb, Tb);
            j = j+1;
        end
        i = i+1;
    end
    j = 1;
    for Re = all_Re
        NRVnoise(j) = solve_NRV(tauf, taud, U, Rn, Re, N, N, Tb)-NRVbasal;
        j = j+1;
    end

    NRV = NRV - NRVbasal;
    
    if syntype == 1
        NRVe_basal = NRVbasal;
        NRVe_abs = NRV;        
    elseif syntype == 2
        NRVi_basal = NRVbasal;
        NRVi_abs = NRV;
    end

    
    %Finds Nb that maximizes NRV
    for j = 1:nCols
        [val,ind] = max( NRV(:,j) );
        nB_max(j) = ind;

        NRV(:,j) = 100*( NRV(:,j) - NRVnoise(j)  )/NRVnoise(j) ;
    end
    
    if syntype == 1
        Gs1 = NRV;  
        nB_maxe = nB_max;
    elseif syntype == 2  
        Gs2 = NRV;
        nB_maxi = nB_max;
    end

    if syntype == 1
        mine = min(min(NRV));
        maxe = max(max(NRV));
    elseif syntype == 2  
        mini = min(min(NRV));
        maxi = max(max(NRV));
    end  
    

    if syntype == 1
        sb(syntype) = subplot('Position',[.07 .79 .2 .155]); 
        imagesc(NRV); hold on;
        set(gca,'YDir','normal'); hold on;
        plot(nB_max,'k','LineWidth',2); 
        colormap(mymap); caxis([-90 120]);
        h = gca();
        h.XTick = linspace(1,nCols,5);
        h.YTick = linspace(1,nRows,5);
        h.XTickLabel = strread(num2str([2 4 6 8 10]),'%s');
        h.YTickLabel = cellstr(num2str(round(log10([1,10,100,1000,10000]')), '10^%d'));
        xlabel('r_{\delta} [% of r_{bas}]');
        ylabel('N_{ext}');
        title(syntitle);
        
        
        sb(4) = subplot('Position',[.08 .28 .2 .15]);
        plot(all_Re/800, all_Nb(nB_max),'b','LineWidth',2); hold on;
        opt_fe = mean( all_Re./all_Nb(nB_max) );
        
    elseif syntype == 2
        sb(syntype) = subplot('Position',[.07 .535 .2 .155]); 
        imagesc(NRV); hold on;
        set(gca,'YDir','normal'); hold on;
        plot(nB_max,'k','LineWidth',2); 
        colormap(mymap); caxis([-90 120]);
        h = gca();
        h.XTick = linspace(1,nCols,5);
        h.YTick = linspace(1,nRows,5);
        h.XTickLabel = strread(num2str([2 4 6 8 10]),'%s');
        h.YTickLabel = cellstr(num2str(round(log10([1,10,100,1000,10000]')), '10^%d'));
        xlabel('r_{\delta} [% of r_{bas}]');
        ylabel('N_{ext}');
        title(syntitle);
        %text(-.35,1.1,'B','FontSize',13,'FontName','Arial','Fontweight','Bold','Units','normalized');
        
        %plot(sb(4),all_Re, all_Nb(nB_max),'r','LineWidth',2); 
        opt_fi = mean( all_Re./all_Nb(nB_max) );
        
    end
    

end


%Finds Nb that maximizes combined FFI NRV
NRVc = Gs1 - Gs2;

for j = 1:nCols
    [val,ind] = max( NRVc(:,j) );
    nB_max(j) = ind;
    
    optval_ce(j) = Gs1(ind,j);
    optabs_ce(j) = NRVe_abs(ind,j);
    optabs_ee(j) = NRVe_abs(nB_maxe(j),j);
    scatter_ee(j) = NRVe_abs(end,j);
    
    optval_ci(j) = Gs2(ind,j);
    optabs_ci(j) = NRVi_abs(ind,j);
    optabs_ii(j) = NRVi_abs(nB_maxi(j),j);
    scatter_ii(j) = NRVi_abs(end,j);
    
    optval_cc(j) = NRVc(ind,j);
    
    %NRVc(:,j) = NRVc(:,j) - NRVc(end,j);
end

% minc = min(min(NRVc));
% maxc = max(max(NRVc));
% mincolor = min([mini, mine, minc]);
% maxcolor = max([maxi, maxe, maxc]);
% set(sb(1),'clim',[mincolor,maxcolor]);
% set(sb(2),'clim',[mincolor,maxcolor]);

sb(3) = subplot('Position',[.35 .535 .25 .41]); 
imagesc(NRVc);
set(gca,'YDir','normal'); hold on;
plot(nB_max,'k','LineWidth',2);
colormap(mymap); caxis([-90 120]); c = colorbar();
set(c,'Position',[.605 .535 .02 .41]);
%set(sb(3),'clim',[mincolor,maxcolor]);
h = gca();
h.XTick = linspace(1,nCols,5);
h.YTick = linspace(1,nRows,5);
h.XTickLabel = strread(num2str([2 4 6 8 10]),'%s');
h.YTickLabel = cellstr(num2str(round(log10([1,10,100,1000,10000]')), '10^%d'));
xlabel('r_{\delta} [% of r_{bas}]');
ylabel('N_{ext}');
title('combined');
title(c,'G^{com}[%]');


plot(sb(1),nB_max,'--k','LineWidth',2);
plot(sb(2),nB_max,'--k','LineWidth',2);

axes(sb(4));
all_Re = all_Re/800;
plot(all_Re, all_Nb(nB_max),'k','LineWidth',2); box off;
xlim([min(all_Re),max(all_Re)]);
h.XTick = linspace(1,nCols,5);
h.XTickLabel = strread(num2str([2 4 6 8 10]),'%s');
opt_fc = mean( all_Re*800./all_Nb(nB_max) );
ylabel('N_{opt}');
%xlabel('r_{\delta} [% of r_{bas}]');
l1 = legend({['r_{opt}^{s1} = ',num2str(round(opt_fe)),' Hz'],...
        ['r_{opt}^{com} = ',num2str(round(opt_fc)),' Hz']});
legend boxoff;



sb(5) = subplot('Position',[.08 .08 .2 .15]);
plot(all_Re,optval_cc,'k','LineWidth',2); hold on;
plot(all_Re,optval_ce,'--b','LineWidth',2);
plot(all_Re,optval_ci,'--r','LineWidth',2); 
plot(all_Re,zeros(1,length(all_Re)),'color',[.5, .5, .5],'LineWidth',.7); box off;
xlim([min(all_Re),max(all_Re)]);
ylabel('G_{max}^{com}[%]');
xlabel('r_{\delta} [% of r_{bas}]');



sb(6) = subplot('Position',[.72 .55 .24 .35]); 
h1 = semilogx(all_Nb,Gs1(:,121),'b','LineWidth',2); hold on;
h2 = semilogx(all_Nb,Gs2(:,121),'r','LineWidth',2);
h3 = semilogx(all_Nb,Gs1(:,121)-Gs2(:,121),'color',[0,0,0],'LineWidth',2);
semilogx(all_Nb,zeros(1,length(all_Nb)),'--k'); box off;
[val, indy] = max( Gs1(:,121)-Gs2(:,121) );
hmax = plot(all_Nb(indy), val, 'ok','LineWidth',2);
text(all_Nb(indy)+20, val+10,['N_{opt}^{com} = ',num2str(round(all_Nb(indy)))]);
l2 = legend([h1 h2 h3],'G^{s1}','G^{s2}','G^{com}'); legend boxoff;
xlabel('N_{ext}'); ylabel('G [%]');


plot(sb(1),121*ones(1,nRows),1:nRows,'color',[.6, .6, .6],'LineWidth',2);
plot(sb(2),121*ones(1,nRows),1:nRows,'color',[.6, .6, .6],'LineWidth',2);
plot(sb(3),121*ones(1,nRows),1:nRows,'color',[.6, .6, .6],'LineWidth',2);


%--------------------------------------------------------------------------
%Second figure for joint optimization (s1 - s2)
%--------------------------------------------------------------------------
gendata = 0;
if gendata == 1
    Rn = .5; 
    Tb = 0.04;
    kEI = 1;

    nSyn = 200;
    all_U = linspace(.05,.9,nSyn);
    all_taud = linspace(.02,.5,nSyn);
    all_tauf = fliplr(all_taud);

    Rs = zeros(nSyn,nSyn);
    Gn = zeros(nSyn,nSyn);   
    for syn1 = 1:nSyn
        U1 = all_U(syn1);
        tauf1 = all_tauf(syn1);
        taud1 = all_taud(syn1);
        for syn2 = 1:nSyn
            U2 = all_U(syn2);
            tauf2 = all_tauf(syn2);
            taud2 = all_taud(syn2);

%             [optfreq, optgain] = theoretical_optfreq_joint(tauf1, taud1, U1, tauf2, taud2, U2, kEI, Rn, Tb);  
%             Rs(syn1,syn2) = optfreq;
%             Gn(syn1,syn2) = optgain;
            
            [Freq, Gain] = search_optfreq_joint(tauf1, taud1, U1, tauf2, taud2, U2, kEI, Rn, Tb);
            Rs(syn1,syn2) = Freq;
            Gn(syn1,syn2) = Gain;
        end
    end
    %if the gain curves are too similar (e.g. syn1=syn2 and kei=1), the combined gain curve doesn't
    %make sense anymore and any Rs value will give the same gain=0.
    Rs(Rs==max(max(Rs))) = 0;
    %save('data_fig5_1kei.mat','Rs','Gn');
end


load('data_fig5_1kei.mat');

sb(7) = subplot('Position',[.37 .07 .255 .35]); 
imagesc(Gn); hold on;
plot(154, 13, 'ok','LineWidth',2);
%plot(,'s','MarkerSize',9,'Linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[1,1,1]);
set(gca,'YDir','normal');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
str1 = '$$ s1: Facilitation \longleftrightarrow  Depression $$';
str2 = '$$ s2: Facilitation \longleftrightarrow  Depression $$';
ylabel(str1,'Interpreter','latex');
xlabel(str2,'Interpreter','latex');
%colormap(gca,'jet'); 
colormap(mymap); caxis([-90 120]); 
c2 = colorbar(); set(c2, 'ylim', [0 120])
set(c2,'Position',[.63 .07 .02 .35]);
title(c2,'G_{max}^{com}[%]');


%--------------------------------------------------------------------------
% Effects of Ongoing basal rate
%--------------------------------------------------------------------------
kEI = 1;
tauf1 = .2;
taud1 = .05;
all_U1 = [0.05, 0.075, 0.1];
tauf2 = .05;
taud2 = .2;
U2 = .7;
all_noise = linspace(0.1,10,60);

Freq = zeros(1,length(all_noise));
Gain = zeros(1,length(all_noise));
i = 0;
for U1 = all_U1
    i = i + 1;
    j = 0;
    for Rn = all_noise  
        j = j+1;
        
        [optfreq, optgain] = theoretical_optfreq_joint(tauf1, taud1, U1, tauf2, taud2, U2, kEI, Rn, Tb);

        Freq(i,j) = optfreq;
        Gain(i,j) = optgain;
    end
end

sb(8) = subplot('Position',[0.65 0.45 0.22 0.10]);
plot(all_noise,Freq(1,:),'Linewidth',2,'color',[.8,.8,.8]); hold on; 
plot(all_noise,Freq(2,:),'Linewidth',2,'color',[.4,.4,.4]);
plot(all_noise,Freq(3,:),'Linewidth',2,'color',[0,0,0]); box off;
ylim([70, 230]);
ylabel('r_{opt}^{com} [Hz]'); set(gca,'xticklabel',[]);
l3 = legend('U^{s1}=0.050','U^{s1}=0.075','U^{s1}=0.100');
legend boxoff;

sb(9) = subplot('Position',[0.65 0.28 0.22 0.10]);
plot(all_noise,Gain(1,:),'Linewidth',2,'color',[.8,.8,.8]); hold on; 
plot(all_noise,Gain(2,:),'Linewidth',2,'color',[.4,.4,.4]);
plot(all_noise,Gain(3,:),'Linewidth',2,'color',[0,0,0]); box off;
ylim([0, 250]);
ylabel('G_{max}^{com} [%]'); set(gca,'xticklabel',[]);


optimum_combined;

sb(10) = subplot('Position',[0.65 0.08 0.22 0.10]);
plot(all_noise,Gain2(1,:),'b','Linewidth',2); hold on; 
plot(all_noise,Gain2(2,:),'r','Linewidth',2);
plot(all_noise,Gain2(3,:),'k','Linewidth',2);
plot([0, max(all_noise)],[0, 0],'k'); box off;
plot([.5,.5],[-100,200],'--k');
xlim([0,max(all_noise)]); ylim([-100, 160]);
xlabel('r_{bas} [Hz]'); ylabel('G_{max} [%]');
l4 = legend('s1','s2','com'); legend boxoff;


%--------------------------------------------------------------------------
% General figure configuration
%--------------------------------------------------------------------------
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',12);
set(findall(fig,'-property','FontName'),'FontName','Arial');
set(l1,'Fontsize',10,'Position',[0.09 0.37 0.13 0.06]);
set(l2,'Fontsize',10,'Position',[0.86 0.83 0.16 0.06]);
set(l3,'Fontsize',10,'Position',[0.85 0.515 0.13 0.06]);
set(l4,'Fontsize',10,'Position',[0.895 0.175 0.08 0.06]);


% text(sb(1),-.32,1.1,'A','FontSize',13,'FontName','Arial','Fontweight','Bold','Units','normalized');
% text(sb(3),-.22,1.05,'B','FontSize',13,'FontName','Arial','Fontweight','Bold','Units','normalized');
% text(sb(4),-.3,1.1,'C','FontSize',13,'FontName','Arial','Fontweight','Bold','Units','normalized');
% text(sb(5),-.3,1.1,'D','FontSize',13,'FontName','Arial','Fontweight','Bold','Units','normalized');
% text(sb(6),-.15,1.08,'E','FontSize',12,'FontName','Arial','Fontweight','Bold','Units','normalized');
% text(sb(7),-.1,1.08,'F','FontSize',12,'FontName','Arial','Fontweight','Bold','Units','normalized');

set(sb(1),'Position',[.07 .79 .2 .155]); %s1 surface
set(sb(2),'Position',[.07 .535 .2 .155]); %s2 surface
set(sb(3),'Position',[.35 .535 .25 .41]); %s1-s2 surface
set(sb(4),'Position',[.08 .28 .2 .15]); %Nopt vs Rext
set(sb(5),'Position',[.08 .08 .2 .15]); %Gopt vs Rext
set(sb(6),'Position',[.73 .69 .23 .25]); %Gain curves
set(sb(7),'Position',[.34 .07 .255 .35]); %Gain surface fac->dep
set(c2,'Position',[.60 .07 .02 .35]);
set(sb(8),'Position',[0.745 0.425 0.215 0.145]);
set(sb(9),'Position',[0.745 0.255 0.215 0.145]);
set(sb(10),'Position',[0.745 0.08 0.215 0.145]);
