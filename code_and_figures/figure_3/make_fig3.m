% Effects of different Ts on Ropt and Gopt
clear all; %close all;

fig1 = figure(); set(gcf,'color','w','Position', [50, 50, 1000, 700]);

%Rn = .5; 
%Tb = 0.3;

nSyn = 60;
all_U = linspace(.05,.9,nSyn);
all_taud = linspace(.02,.5,nSyn);
all_tauf = fliplr(all_taud);

Rs = zeros(1,nSyn);
PRRp = zeros(1,nSyn);
all_Rs = [];
all_PRRp = [];
i = 0;
for Tb = [0.04, 0.1, 0.3] 
    for Rn = [.5, 2, 5]
        i = i + 1;
        for syntype = 1:length(all_U)
            U = all_U(syntype);
            tauf = all_tauf(syntype);
            taud = all_taud(syntype);

            [optfreq, optgain] = theoretical_optfreq(tauf, taud, U, Rn, Tb);

            Rs(syntype) = optfreq;
            PRRp(syntype) = mean( optgain );

        end
        
        all_Rs(i,:) = Rs;
        all_PRRp(i,:) = PRRp;

        if Tb == 0.04
            if  Rn == .5%Tb == 0.04
                sb1 = subplot('Position',[0.09 0.73 0.38 0.22]); 
                plot(Rs,'color',[0,0,0],'LineWidth',2); hold on;
                sb2 = subplot('Position',[0.09 0.45 0.38 0.22]);
                plot(PRRp,'color',[0,0,0],'LineWidth',2); hold on;
                
            elseif Rn == 2%Tb == 0.1
                plot(sb1,Rs,'color',[.5,.5,.5],'LineWidth',2);
                plot(sb2,PRRp,'color',[.5,.5,.5],'LineWidth',2);

            elseif Rn == 5%Tb == 0.3
                plot(sb1,Rs,'color',[.8,.8,.8],'LineWidth',2);
                plot(sb2,PRRp,'color',[.8,.8,.8],'LineWidth',2);

            end
        end
    end
end


set(fig1, 'CurrentAxes', sb1);
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
xlim([0,nSyn+1]); ylim([0, 200]); box off;
plot(5,0,'ob','Linewidth',3,'Markersize',10);
plot(47,0,'or','Linewidth',3,'Markersize',10);
l4 = legend({'r_{bas} = 0.5 Hz','r_{bas} = 2 Hz','r_{bas} = 5 Hz','s1','s2'}); 
legend boxoff;
%set(l4,'Position',[0.85 0.35 0.13 0.11]);
ylabel('r_{opt} [Hz]');
str1 = '$$ Facilitation \longleftrightarrow  Depression $$';
xlabel(str1,'Interpreter','latex');



set(fig1, 'CurrentAxes', sb2);
set(sb2, 'XTickLabelMode', 'manual', 'XTickLabel', []);
%set(gca,'Position',[0.78 0.1 0.19 0.16]);
ylabel('G_{max} [%]');
str1 = '$$ Facilitation \longleftrightarrow  Depression $$';
xlabel(str1,'Interpreter','latex');
xlim([0,nSyn+1]); ylim([0, 150]); box off;
plot(5,0,'ob','Linewidth',3,'Markersize',10);
plot(47,0,'or','Linewidth',3,'Markersize',10);


sb3 = subplot('Position',[0.09 0.09 0.38 0.27]); 
plot(all_Rs(1,:),all_PRRp(1,:),'Linewidth',2,'color',[0,0,0]); hold on;
plot(all_Rs(2,:),all_PRRp(2,:),'Linewidth',2,'color',[.5,.5,.5]);
plot(all_Rs(3,:),all_PRRp(3,:),'Linewidth',2,'color',[.8,.8,.8]); box off;
xlim([0,200]); ylim([0,150]); 
xlabel('r_{opt} [Hz]');
ylabel('G_{max} [%]');
str2 = '$$ Depression \longleftrightarrow  Facilitation $$';
text('String',str2,'Position',[46 30 0],'Rotation',25,'Interpreter','latex','FontName','Arial');





rd = 0.1;
sb4 = subplot('Position',[0.57 0.73 0.38 0.22]); 
plot(rd./all_Rs(9,:),'-.','Linewidth',2,'color',[0.8, 0.8, 0.8]); hold on;
plot(rd./all_Rs(6,:),'--','Linewidth',2,'color',[0.8, 0.8, 0.8]); 
plot(rd./all_Rs(3,:),'Linewidth',2,'color',[0.8, 0.8, 0.8]);
plot(rd./all_Rs(8,:),'-.','Linewidth',2,'color',[0.5, 0.5, 0.5]);
plot(rd./all_Rs(5,:),'--','Linewidth',2,'color',[0.5, 0.5, 0.5]);
plot(rd./all_Rs(2,:),'Linewidth',2,'color',[0.5, 0.5, 0.5]);
pl3 = plot(rd./all_Rs(7,:),'-.','Linewidth',2,'color',[0, 0, 0]); 
pl2 = plot(rd./all_Rs(4,:),'--','Linewidth',2,'color',[0, 0, 0]); 
pl1 = plot(rd./all_Rs(1,:),'Linewidth',2,'color',[0, 0, 0]); box off;
l5 = legend([pl1,pl2,pl3],{'T_s = 40 ms','T_s = 100 ms','T_s = 300 ms'},'AutoUpdate','off'); 
legend boxoff;
plot(5,0,'ob','Linewidth',3,'Markersize',10);
plot(47,0,'or','Linewidth',3,'Markersize',10);
xlim([0,nSyn+1]);
ylabel('OD');
xlabel(str1,'Interpreter','latex');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);




load('all_STP_40ms.mat','Rs','PRRp','sparseness');
nSyn = 20;
all_U = linspace(.05,.9,nSyn);
all_taud = linspace(.02,.5,nSyn);
all_tauf = fliplr(all_taud);

X = repelem(all_taud,nSyn*nSyn);
Y = repmat(repelem(all_tauf,nSyn),1,nSyn);
Z = repmat(all_U,1,nSyn*nSyn);
S = 20./sparseness;

sb5 = subplot('Position',[0.57 0.11 0.38 0.5]); 
scatter3(X,Y,Z,S,'filled');
xlabel('\tau_{rec}');
ylabel('\tau_{f}');
zlabel('U');
xlim([0,0.5]); ylim([0,0.5]); zlim([0,0.9]);
campos([3.25 -2.75 2]);
l6 = legend('OD'); legend boxoff;
set(l6,'Position',[0.83 0.6 0.15 0.04]);



fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',12);
set(findall(fig,'-property','FontName'),'FontName','Arial');