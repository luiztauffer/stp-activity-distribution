%Generate figure
clear all; close all;

load('fig4data_40ms.mat');

all_taud = linspace(0.01,0.1,600);
all_U = linspace(0.05,0.15,161);

nRows = length(all_taud);
nCols = length(all_U);

for i = 1:nCols
    [val,ind] = min( abs(Gain(:,i)-108.4) );
    if val < 1
        G75x(i) = i;
        G75y(i) = ind;
    end
    
    
    [val,ind] = min( abs(Freq(:,i)-150) );
    if val < 1
        F150x(i) = i;
        F150y(i) = ind;
        GF150(i) = Gain(ind,i);
        TD2(1,i) = all_taud(ind);
    end
    
    [val,ind] = min( abs(Freq(:,i)-175) );
    if val < 1
        F180x(i) = i;
        F180y(i) = ind;
        GF180(i) = Gain(ind,i);
        TD2(2,i) = all_taud(ind);
    end
    
    [val,ind] = min( abs(Freq(:,i)-200) );
    if val < 1
        F200x(i) = i;
        F200y(i) = ind;
        GF200(i) = Gain(ind,i);
        TD2(3,i) = all_taud(ind);
    end
end

for f = 150:200
    [val1,ind1] = min( abs(Freq(:,1)-f) );
    [val2,ind2] = min( abs(Freq(:,41)-f) );
    [val3,ind3] = min( abs(Freq(:,81)-f) );
    
    GF(1,f-149) = Gain(ind1,1);
    GF(2,f-149) = Gain(ind2,41);
    GF(3,f-149) = Gain(ind3,81);
    
    TD1(1,f-149) = all_taud(ind1);
    TD1(2,f-149) = all_taud(ind2);
    TD1(3,f-149) = all_taud(ind3);
end

TD1max = max(max(TD1));
TD1scaled = TD1/TD1max;
TD2max = max(max(TD2));
TD2scaled = TD2/TD2max;

figure(); set(gcf,'color','w','Position', [50, 50, 900, 500]);
sb1 = subplot('Position',[.08 .59 .4 .34]);
imagesc(Freq); hold on;
set(gca,'YDir','normal');
colormap('jet'); c1 = colorbar();
plot(F150x,F150y,'k','Linewidth',2);
plot(1,F150y(1),'o','MarkerSize',9,'Linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[1,1,1]);
plot(81,F150y(81),'s','MarkerSize',9,'Linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[1,1,1]);
title(c1,'r_{opt}[Hz]');
set(c1,'Position',[.49 .59 .02 .34]);
h = gca();
h.XTick = linspace(1,nCols,3);
h.YTick = linspace(1,nRows,4);
h.XTickLabel = strread(num2str(all_U(round(linspace(1,nCols,3)))),'%s');
h.YTickLabel = strread(num2str(round(1000*all_taud(round(linspace(1,nRows,4))))),'%s');
xlabel('U');
ylabel('\tau_{rec} [ms]');


sb2 = subplot('Position',[.08 .1 .4 .34]);
imagesc(Gain); hold on;
set(gca,'YDir','normal');
colormap('jet'); c2 = colorbar();
plot(G75x,G75y,'k','Linewidth',2);
plot(1,F150y(1),'o','MarkerSize',9,'Linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[1,1,1]);
plot(81,F150y(81),'s','MarkerSize',9,'Linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[1,1,1]);
title(c2,'G_{max}[%]');
set(c2,'Position',[.49 .1 .02 .34]);
h = gca();
h.XTick = linspace(1,nCols,3);
h.YTick = linspace(1,nRows,4);
h.XTickLabel = strread(num2str(all_U(round(linspace(1,nCols,3)))),'%s');
h.YTickLabel = strread(num2str(round(1000*all_taud(round(linspace(1,nRows,4))))),'%s');
xlabel('U');
ylabel('\tau_{rec} [ms]');


sb3 = subplot('Position',[.64 .59 .3 .34]);
plot(150:200,GF(1,:),'Linewidth',1,'Color',[.8,.8,.8]); hold on;
plot(150:200,GF(2,:),'Linewidth',1,'Color',[.5,.5,.5]); 
plot(150:200,GF(3,:),'Linewidth',1,'Color',[0,0,0]);
for i = 1:length(GF(1,:))
    plot(149+i,GF(1,i),'o','MarkerSize',10*TD1scaled(1,i),'Color',[.8,.8,.8],'MarkerFaceColor',[.8,.8,.8]); hold on;
    plot(149+i,GF(2,i),'o','MarkerSize',10*TD1scaled(1,i),'Color',[.5,.5,.5],'MarkerFaceColor',[.5,.5,.5]); 
    plot(149+i,GF(3,i),'o','MarkerSize',10*TD1scaled(1,i),'Color',[0,0,0],'MarkerFaceColor',[0,0,0]); 
end
box off;
xlim([145, 205]); ylim([80,150]);
xlabel('r_{opt}[Hz]');
ylabel('G_{max}[%]');
l1 = legend('U = 0.050','U = 0.075','U = 0.100'); legend boxoff;
set(l1,'Position',[0.66 0.84 0.13 0.12]);


sb4 = subplot('Position',[.64 .1 .3 .34]);
plot(all_U(F150x),GF150,'Linewidth',1,'color',[.8,.8,.8]); hold on; 
plot(all_U(F180x),GF180,'Linewidth',1,'color',[.5,.5,.5]);
plot(all_U(F200x),GF200,'Linewidth',1,'color',[0,0,0]);
for i = 1:length(GF150)
    plot(all_U(i),GF150(i),'o','MarkerSize',10*TD2scaled(1,i),'Color',[.8,.8,.8],'MarkerFaceColor',[.8,.8,.8]); hold on; 
end
for i = 1:length(GF180)
    plot(all_U(i),GF180(i),'o','MarkerSize',10*TD2scaled(2,i),'Color',[.5,.5,.5],'MarkerFaceColor',[.5,.5,.5]);
end
for i = 1:length(GF200)    
    plot(all_U(i),GF200(i),'o','MarkerSize',10*TD2scaled(3,i),'Color',[0,0,0],'MarkerFaceColor',[0,0,0]); 
end
box off;
l2 = legend('150 Hz','175 Hz','200 Hz'); legend boxoff;
%set(l2,'Position',[0.66 0.84 0.13 0.12]);
ylim([80, 150]); %xlim([0.045, 0.105]);
xlabel('U');
ylabel('G_{max}[%]');


fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',12);
set(findall(fig,'-property','FontName'),'FontName','Arial');

% text(sb1,-0.17,1.05,'A','FontSize',12,'FontName','Arial','Fontweight','Bold','Units','normalized');
% text(sb2,-0.17,1.05,'B','FontSize',12,'FontName','Arial','Fontweight','Bold','Units','normalized');
% text(sb3,-0.19,1.05,'C','FontSize',12,'FontName','Arial','Fontweight','Bold','Units','normalized');
% text(sb4,-0.19,1.05,'D','FontSize',12,'FontName','Arial','Fontweight','Bold','Units','normalized');

%print('fig41','-dsvg');



