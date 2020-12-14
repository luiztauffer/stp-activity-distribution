clear all; close all;

% Shematic figure
shift_1 = 20;slope_1 = 0.04;
syn_amp = [0:1:160];
fr_sch_1 = .8./(1+(shift_1*exp(-slope_1*(syn_amp+5))));

shift_21 = 20;slope_21 = 0.08;
fr_sch_21 = 0.3./(1+(shift_21*exp(-slope_21*(syn_amp-40))));
shift_22 = 10;slope_22 = 0.3;
fr_sch_22 = 0.2./(1+(shift_22*exp(-slope_22*(syn_amp-60))));

nmda_tx = (fr_sch_22+fr_sch_21+fr_sch_1*0.7)+0.04;
nonmda_tx = fr_sch_1;
max_val = max(nmda_tx);
nmda_tx1 = nmda_tx/max_val;
nonmda_tx = nonmda_tx/max_val;
nmda_tx = [nonmda_tx(1:40) nmda_tx1(41:end)];

fig = figure(); 
set(gcf,'color','w','Position', [50, 50, 1000, 400]);
% clf
dark_Nonmda = [0.1 0.1 0.6];
dark_nmda = [0.6 0.1 0.1];
light_Nonmda = [0.7 0.7 0.9];
light_nmda = [0.9 0.7 0.7];
Nonmda_col = [0.2 0.2 0.8];
nmda_col = [0.8 0.2 0.2];
% 

p1 = plot(syn_amp,nmda_tx,'linewidth',2,'DisplayName','With NMDA','color',nmda_col);
hold on
p2 = plot(syn_amp,nonmda_tx,'linewidth',2,'DisplayName','With NMDA','color',Nonmda_col);

input_val = [90 130];

spk_out_nmda = interp1(syn_amp,nmda_tx,input_val);
spk_out_Nonmda = interp1(syn_amp,nonmda_tx,input_val);

p3 = plot([input_val(1) input_val(1)],[0 spk_out_nmda(1)],'--','linewidth',1,'color',dark_nmda);
p4 = plot([input_val(2) input_val(2)],[0 spk_out_nmda(2)],'--','linewidth',1,'color',dark_nmda);
p5 = plot([0 input_val(1)], [spk_out_nmda(1) spk_out_nmda(1)],'--','linewidth',1,'color',dark_nmda);
p6 = plot([0 input_val(2)], [spk_out_nmda(2) spk_out_nmda(2)],'--','linewidth',1,'color',dark_nmda);

p7 = plot([0 input_val(1)], [spk_out_Nonmda(1) spk_out_Nonmda(1)],'--','linewidth',1,'color',dark_Nonmda);
p8 = plot([0 input_val(2)], [spk_out_Nonmda(2) spk_out_Nonmda(2)],'--','linewidth',1,'color',dark_Nonmda);

py1 = plot(input_val(1),spk_out_nmda(1),'.','markersize',20,'color',dark_nmda);
py2 = plot(input_val(2),spk_out_nmda(2),'.','markersize',20,'color',dark_nmda);
py3 = plot(input_val(1),spk_out_Nonmda(1),'.','markersize',20,'color',dark_Nonmda);
py4 = plot(input_val(2),spk_out_Nonmda(2),'.','markersize',20,'color',dark_Nonmda);

py1 = plot(input_val(1),spk_out_nmda(1),'.','markersize',10,'color',[0.8 0.8 0.8]);
py2 = plot(input_val(2),spk_out_nmda(2),'.','markersize',10,'color',[0.8 0.8 0.8]);
py3 = plot(input_val(1),spk_out_Nonmda(1),'.','markersize',10,'color',[0.8 0.8 0.8]);
py4 = plot(input_val(2),spk_out_Nonmda(2),'.','markersize',10,'color',[0.8 0.8 0.8]);

toh = text(1,1,'3');
set(toh,'unit','data','position',[input_val(1)+2 spk_out_nmda(1)-.05],'color',dark_nmda,'fontsize',20,'fontweight','bold');
toh = text(1,1,'4');
set(toh,'unit','data','position',[input_val(2)+2 spk_out_nmda(2)-.05],'color',dark_nmda,'fontsize',20,'fontweight','bold');
toh = text(1,1,'1');
set(toh,'unit','data','position',[input_val(1)+2 spk_out_Nonmda(1)-.05],'color',dark_Nonmda,'fontsize',20,'fontweight','bold');
toh = text(1,1,'2');
set(toh,'unit','data','position',[input_val(2)+2 spk_out_Nonmda(2)-.05],'color',dark_Nonmda,'fontsize',20,'fontweight','bold');

set(gca,'XTick',[0:40:160],'XTickLabel',syn_amp(1:40:end)/160);
set(gca,'YTick',[0:0.2:1]); %'XTickLabel',syn_amp(1:40:end)/160);

% legend([p1 p2],{'With NMDA','Without NMDA'},'Location','NorthWest','NumColumns',2)
legend([p2 p1], {'TF-1','TF-2'}, 'Location', 'NorthEast');
% legend('boxoff');

y_lab_txt = sprintf('%s','Output spike prob.');
% set(gca,'Linewidth',2,'fontweight','bold','fontsize',12)
xlabel('Efefctive input == Q^{S1} (a.u.)')
ylabel(y_lab_txt);
ylim([0 1.1]);
xlim([0 160]);

text(0.58*160, 0.1, 'dense');
text(0.82*160, 0.1, 'sparse');

set(gca,'Position',[.1 .16 .7 .8]);


% set(gcf, 'PaperUnits','centimeters');
% set(gcf, 'PaperSize',[12 12]);
% set(fig,'PaperPosition',[0,0,(get(fig,'PaperSize'))])
box off;
set(findall(fig,'-property','FontSize'),'FontSize',12);
set(findall(fig,'-property','FontName'),'FontName','Arial');



