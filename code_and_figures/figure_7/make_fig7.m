clear all; close all;

load('Expected_Gain.mat')
nShape = length(GammaShape);
nScale = length(GammaScale);

Gain1 = interp2(Gain1, 'linear');
Gain2 = interp2(Gain2, 'linear');
xq = 1.5:(nShape+.5);
GammaShape_in = interp1(GammaShape, xq, 'linear');
GammaShape_extended = reshape([GammaShape;GammaShape_in],1,[]);
GammaShape_extended(end) = [];
xq = 1.5:(nScale+.5);
GammaScale_in = interp1(GammaScale, xq, 'linear');
GammaScale_extended = reshape([GammaScale;GammaScale_in],1,[]);
GammaScale_extended(end) = [];

nShape = length(GammaShape_extended);
nScale = length(GammaScale_extended);

nRates = 100000;
rates = linspace(0,500,nRates);
optDist = zeros(nShape,nRates);

%Combined Gain
Gainc = Gain1-Gain2;

for i = 1:nShape
    %Location of max Gain per row
    [~, max_loc1(i)] = max(Gainc(i,:));
    [~, max_loc2(i)] = max(Gainc(i,:));
    
    optShape(i) = GammaShape_extended(i);
    optScale(i) = GammaScale_extended(max_loc1(i));
    
    uRopt(i) = optShape(i)*optScale(i);
    
    optDist(i,:) = gampdf(rates,optShape(i),optScale(i));
    [~, uropt_ind(i)] = min(abs(uRopt(i)-rates));   %mean ropt for Gamma dist
    [~, biropt_ind(i)] = min(abs(131-rates));  %ropt for bimodal dist
end

sh_ind = [1,ceil(nShape/2),nShape];
sc_ind = [nScale-100,max_loc1(1),100;...
          nScale-100,max_loc1(ceil(nShape/2)),100;...
          nScale-100,max_loc1(nShape),100];

mymap = my_cmap(-90,120,0);   %Blue to Red, with center on 0

figure(); set(gcf,'color','w','Position', [50, 50, 1000, 330]);
sb1 = subplot(1,3,1);
imagesc(Gain1); hold on;
colormap(mymap); caxis([-90 120]);
set(gca,'XDir','reverse');
plot([1,nScale],[ceil(nShape/2),ceil(nShape/2)],'color',[.7,.7,.7],'Linewidth',1.5);
plot(fliplr(max_loc1),linspace(nShape,1,nShape),'k','Linewidth',1.5);
h = gca();
h.YTick = round(linspace(1,nShape,5));
h.XTick = round(linspace(1,nScale,6));
h.YTickLabel = strread(num2str([1 5 10 15 20]),'%s');
h.XTickLabel = cellstr(num2str(round(log10([.01,.1,1,10,100,1000]')),'10^{%1.f}'));
ylabel('Shape','FontSize',12); xlabel('Scale','FontSize',12);
title('s1 - facilitatory','FontSize',12);

sb2 = subplot(1,3,2);
imagesc(Gain2); hold on;
plot([1,nScale],[ceil(nShape/2),ceil(nShape/2)],'color',[.7,.7,.7],'Linewidth',1.5);
plot(fliplr(max_loc2),linspace(nShape,1,nShape),'k','Linewidth',1.5);
colormap(mymap); caxis([-90 120]);
%colormap(brewermap([],'*RdBu')); c2 = colorbar(); caxis([-90 60]);
set(gca,'XDir','reverse');
h = gca();
h.YTick = round(linspace(1,nShape,5));
h.XTick = round(linspace(1,nScale,6));
h.YTickLabel = strread(num2str([1 5 10 15 20]),'%s');
h.XTickLabel = cellstr(num2str(round(log10([.01,.1,1,10,100,1000]')), '10^{%1.f}'));
ylabel('Shape','FontSize',12); xlabel('Scale','FontSize',12);
title('s2 - depressing','FontSize',12);

sb3 = subplot(1,3,3);
imagesc(Gainc); hold on;
plot([1,nScale],[ceil(nShape/2),ceil(nShape/2)],'color',[.7,.7,.7],'Linewidth',1.5);
plot(fliplr(max_loc2),linspace(nShape,1,nShape),'k','Linewidth',1.5);
colormap(mymap); c2 = colorbar(); caxis([-90 120]);
set(gca,'XDir','reverse');
plot(fliplr(max_loc1),linspace(nShape,1,nShape),'k','Linewidth',1.5);
ii = 0;
for i = sh_ind
    ii = ii+1;
    for j = 1:3
        jj = sc_ind(ii,j);
        if j == 2
            plot(jj,i,'sqk','MarkerSize',11,'Linewidth',4);
        else
            plot(jj,i,'sqk','MarkerSize',11,'Linewidth',1.5);
        end
    end 
end
h = gca();
h.YTick = round(linspace(1,nShape,5));
h.XTick = round(linspace(1,nScale,6));
h.YTickLabel = strread(num2str([1 5 10 15 20]),'%s');
h.XTickLabel = cellstr(num2str(round(log10([.01,.1,1,10,100,1000]')), '10^{%1.f}'));
ylabel('Shape','FontSize',12); xlabel('Scale','FontSize',12);
title('combined','FontSize',12);
title(c2,'G^{com} [%]','FontSize',12);


set(sb1,'Position',[.05 .16 .25 .74]);
set(sb2,'Position',[.365 .16 .25 .74]);
set(sb3,'Position',[.68 .16 .25 .74]);
set(c2,'Position',[.94 .16 .02 .74]);



figure(); set(gcf,'color','w','Position', [50, 50, 630, 410]);
i = 0;
for sh = GammaShape_extended(sh_ind)
    i = i + 1;
    for j = 1:3
        sc = GammaScale_extended(sc_ind(i,j));
        rates_pdf = gampdf(rates,sh,sc);
        
        sb((i-1)*3+j) = subplot(3,3,(i-1)*3+j);
        plot(rates,rates_pdf,'k','Linewidth',1.5);
        if j == 3
            xlim([0,2]);
        else
            xlim([0,400]);
        end
        %Text with the mean rate value
        ylims = get(gca,'ylim');
        xlims = get(gca,'xlim');
        if j==3
            text(xlims(2)*.4,ylims(2),['mean: ',num2str(sh*sc,'%1.2f')],...
                 'VerticalAlignment','top','HorizontalAlignment','left');
        else
            text(xlims(2)*.4,ylims(2),['mean: ',num2str(ceil(sh*sc))],...
                 'VerticalAlignment','top','HorizontalAlignment','left');
        end
        %Y ans X labels
        if (i==2) && (j==1)
            ylabel('P(r_{ext})','FontSize',12);
        elseif (i==3) && (j==2)
            xlabel('r_{ext} [spks/sec]','FontSize',12);
        end
        %Box line width
        if j == 2
            set(gca,'linewidth',3);
        else
            set(gca,'linewidth',.7);
        end
        %Remove xticks
        if i ~= 3
            set(gca,'xtick',[]);
            set(gca,'xticklabel',[]);
        end
    end
end

set(sb(1),'Position',[.1 .72 .23 .21]);
set(sb(2),'Position',[.4 .72 .23 .21]);
set(sb(3),'Position',[.7 .72 .23 .21]);
set(sb(4),'Position',[.1 .44 .23 .21]);
set(sb(5),'Position',[.4 .44 .23 .21]);
set(sb(6),'Position',[.7 .44 .23 .21]);
set(sb(7),'Position',[.1 .16 .23 .21]);
set(sb(8),'Position',[.4 .16 .23 .21]);
set(sb(9),'Position',[.7 .16 .23 .21]);




%Figure 3 -----------------------------------------------------------------
g1curve = Gain1(ceil(nShape/2),:);
g2curve = Gain2(ceil(nShape/2),:);
gccurve = Gainc(ceil(nShape/2),:);

figure(); set(gcf,'color','w','Position', [50, 50, 360, 410]);
h4 = subplot(2,1,1);
plot(g1curve,'b','Linewidth',2); hold on;
plot(g2curve,'r','Linewidth',2);
plot(gccurve,'k','Linewidth',2);
plot([0,nScale],[0,0],'--k'); box off;
h4.XTick = round(linspace(1,nScale,6));
h4.XTickLabel = cellstr(num2str(round(log10([.01,.1,1,10,100,1000]')), '10^{%1.f}'));
set(gca,'XDir','reverse');
ylim([-100, 150]); xlim([0,nScale]); 
ylabel('Gain [%]','Fontsize',12);
xlabel('Scale','Fontsize',12);
l2 = legend({'s1','s2','com'}); legend boxoff;
set(l2,'Position',[0.7 0.8 0.23 0.16],'FontSize',12); legend boxoff;


h5 = subplot(2,1,2);
imagesc(optDist'); hold on; c3 = colorbar(); caxis([0 0.015]);
colormap(brewermap([],'PuBu')); box off;
pline = plot(uropt_ind,'k','Linewidth',2); 
bir = plot(biropt_ind,'--k');
set(gca,'YDir','normal');
h5.YTick = round(linspace(1,50000,5));
h5.XTick = round(linspace(1,nShape,5));
h5.YTickLabel = strread(num2str(round(rates(round(linspace(1,50000,5))))),'%s');
h5.XTickLabel = strread(num2str([1 5 10 15 20]),'%s');
ylim([0,50000]);
ylabel('r_{ext}','Fontsize',12);
xlabel('Shape','Fontsize',12);
title(c3,'P(r_{ext})','FontSize',10);
l1 = legend([pline,bir],{'mean r_{ext}','bi r_{opt}'}); 
set(l1,'Position',[0.52 0.33 0.3 0.10],'FontSize',10); legend boxoff;

set(h4,'Position',[.18 .62 .65 .34]);
set(h5,'Position',[.18 .12 .65 .34]);
set(c3,'Position',[.85 .12 .04 .34]);

