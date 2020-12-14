%Generates figure
clear all; close all;

% Load simulation data
load(['Sim_160000_05/New_Neuron_sim_1_perc_8.mat']);
spkt1 = spkt;
VmP1 = VmP;
load(['Sim_160000_05/New_Neuron_sim_2_perc_8.mat']);
spkt2 = spkt;
VmP2 = VmP;
load(['Sim_160000_05/New_Neuron_sim_3_perc_8.mat']);
spkt3 = spkt;
VmP3 = VmP;
load(['Sim_160000_05/New_Neuron_sim_0_perc_8.mat']);
spkt0 = spkt;
VmP0 = VmP;

nBins = 120;
dt = 0.1;
nTrials = length(VmP)*dt/nBins; %3010;

%spikes within trials ---------------------------------------------
spkt0t = mod(spkt0,120);
spkt1t = mod(spkt1,120);
spkt2t = mod(spkt2,120);
spkt3t = mod(spkt3,120);

edges0 = 0:120;
[spkc0,~] = histcounts(spkt0t,edges0);
[spkc1,~] = histcounts(spkt1t,edges0);
[spkc2,~] = histcounts(spkt2t,edges0);
[spkc3,~] = histcounts(spkt3t,edges0);

%PSTH -------------------------------------------------------------
spkavg0 = 1000*spkc0/nTrials;
spkavg1 = 1000*spkc1/nTrials;
spkavg2 = 1000*spkc2/nTrials;
spkavg3 = 1000*spkc3/nTrials;

%Average rates ----------------------------------------------------
Rbasal = mean(spkavg0(1:40));
Rstim0 = mean(spkavg0(41:80));
Rstim1 = mean(spkavg1(41:80));
Rstim2 = mean(spkavg2(41:80));
Rstim3 = mean(spkavg3(41:80));


%Spk counts per trial ---------------------------------------------
for tr = 1:nTrials
    bas_ini = (tr-1)*120 + eps;
    bas_end = bas_ini + 40;
    stim_ini = bas_end + eps;
    stim_end = stim_ini + 40;

    SpkCnt_b(tr) = sum((spkt0>bas_ini).*(spkt0<bas_end));
    SpkCnt_s0(tr) = sum((spkt0>stim_ini).*(spkt0<stim_end));
    SpkCnt_s1(tr) = sum((spkt1>stim_ini).*(spkt1<stim_end));
    SpkCnt_s2(tr) = sum((spkt2>stim_ini).*(spkt2<stim_end));
    SpkCnt_s3(tr) = sum((spkt3>stim_ini).*(spkt3<stim_end));
end

edges = 0:13;
[cnts_b,~] = histcounts(SpkCnt_b,edges,'Normalization','probability');
[cnts_s0,~] = histcounts(SpkCnt_s0,edges,'Normalization','probability');
[cnts_s3,~] = histcounts(SpkCnt_s3,edges,'Normalization','probability');



%Mean membrane voltage --------------------------------------------
Vm0 = reshape(VmP0,[],nTrials)'; Vm0(1:10,:) = [];
Vm1 = reshape(VmP1,[],nTrials)'; Vm1(1:10,:) = [];
Vm2 = reshape(VmP2,[],nTrials)'; Vm2(1:10,:) = [];
Vm3 = reshape(VmP3,[],nTrials)'; Vm3(1:10,:) = [];
for bin = 1:nBins
    tini = (bin-1)*1/dt + 1;
    tend = tini + 1/dt - 1;

    uVm0(bin) = 1000*mean(mean( Vm0(:,tini:tend) ))-0.2;
    uVm1(bin) = 1000*mean(mean( Vm1(:,tini:tend) ))-0.2;
    uVm2(bin) = 1000*mean(mean( Vm2(:,tini:tend) ))-0.2;
    uVm3(bin) = 1000*mean(mean( Vm3(:,tini:tend) ))-0.2;

    sdVm0(bin) = 1000*std(mean( Vm0(:,tini:tend),2 ));
    sdVm1(bin) = 1000*std(mean( Vm1(:,tini:tend),2 ));
    sdVm2(bin) = 1000*std(mean( Vm2(:,tini:tend),2 ));
    sdVm3(bin) = 1000*std(mean( Vm3(:,tini:tend),2 ));
end



%Plots ------------------------------------------------------------
fig0 = figure(); set(gcf,'color','w','Position', [50, 50, 1200, 300]);

sb1 = subplot(1,3,2);
plot([0, 120],[1,1]*Rbasal,'color',[.3,.3,.3],'Linewidth',.7); hold on;
h11 = plot(spkavg0,'k','Linewidth',2);
%plot(spkavg1);
%plot(spkavg2);
h12 = plot(spkavg3,'color',[.7,.7,.7],'Linewidth',2); box off;
xlim([0,120]); 
xlabel('Time [ms]'); ylabel('Spiking rate [Hz]');
leg1 = legend([h11,h12],{'Sparse';'Smooth'}); legend boxoff;
    

sb2 = subplot(2,3,3);
plot([0, 120],[-53,-53],'color',[.3,.3,.3],'Linewidth',.7); hold on;
h31 = plot(uVm0,'k','Linewidth',2);
h32 = plot(uVm1,'b','Linewidth',2);
h33 = plot(uVm2,'r','Linewidth',2);
h34 = plot(uVm3,'color',[.7,.7,.7],'Linewidth',2); box off;
xlim([0,120]); 
xlabel('Time [ms]'); ylabel('Vm mean [mV]');
leg2 = legend([h31,h32,h33,h34],{'sparse: s1-F, s2-D';'sparse: s1-F, s2-S';'sparse: s1-S, s2-D';'smooth:  s1-F, s2-D'}); legend boxoff;

% sb3 = subplot(2,3,6);
% plot(sdVm0,'k','Linewidth',2); hold on;
% plot(sdVm1,'Linewidth',2);
% plot(sdVm2,'Linewidth',2);
% plot(sdVm3,'color',[.7,.7,.7],'Linewidth',2); box off;
% xlim([0,120]); ylim([0, 2]);
% xlabel('Time [ms]'); ylabel('Vm std [mV]');



%Mean Conductances  --------------------------------------------
fig1 = figure(); set(gcf,'color','w','Position', [50, 50, 1200, 300]);
titulo = {'sparse: s1-F, s2-D';'sparse: s1-F, s2-S';'sparse: s1-S, s2-D';'smooth:  s1-F, s2-D'};
for caso = 1:4
    load(['Sim_160000_05/New_Neuron_sim_',num2str(caso-1),'_perc_8.mat']);
    Ge = reshape(GE,[],nTrials)';
    Gi = reshape(GI,[],nTrials)';

    for bin = 1:nBins
        tini = (bin-1)*1/dt + 1;
        tend = tini + 1/dt - 1;

        uGe(bin) = 1000*mean(mean( Ge(:,tini:tend) ));
        uGi(bin) = 1000*mean(mean( Gi(:,tini:tend) ));
    end

    dGe = 100*(uGe-mean(uGe(1:40)))/mean(uGe(1:40));
    dGi = 100*(uGi-mean(uGi(1:40)))/mean(uGi(1:40));
    
    sbcond(caso) = subplot(1,4,caso);
    plot(dGe,'Linewidth',2); hold on; 
    plot(dGi,'Linewidth',2); 
    plot([1,120],[8,8],'--k'); box off;
    ylim([-.5,16]); xlim([0, 120]);
    yticks([0,4,8,12,16]);
    ylabel('\Delta Conductance [%]'); xlabel('Time [ms]');
    text(30,15.5,titulo{caso});
    
    if caso == 2
        drawbrace([39 8], [39 12.5], 6,'Color','k','Linewidth',1.5);
        text(16,10.8,'G^{s1}');
    elseif caso == 3
        drawbrace([39 2.5], [39 8], 6,'Color','k','Linewidth',1.5);
        text(16,5.8,'G^{s2}');
    end
    
end
lcond = legend({'Exc';'Inh';'Static'}); legend boxoff;



%Plot gradually changing Next ---------------------------------------------
Nopt_perc = [1, 10, 25, 50, 100, 200, 400, 1000, 90000];
load('Varying_Next/PRR_sim_8perc.mat');
all_Nb = all_Nb + 1;
[v,ind] = max(GcC);
Nopt = all_Nb(ind);
i = 0;
for ni = Nopt_perc
    i = i + 1;
    load(['Varying_Next/Nopt_ratio_',num2str(ni),'_perc_8.mat']);

    nTrials = length(VmP)/1200 - 10;
    aux0 = VmP(12001:end);
    aux1 = reshape(aux0,1200,[]);
    aux2 = mean(aux1,2);
    uVm(i) = mean(aux2(401:800));
    
    aux3 = spkt; %round(spkt0*1000);
    aux3(aux3<1200) = [];
    aux4 = mod(aux3,120);
    aux5 = histcounts(aux4,0:10:120);
    uSpks(i) = sum(aux5(5:8))/nTrials;
    
    [~,aux_ind(i)] = min( abs(all_Nb-round(ni*Nopt/100)) );
end
uVm(9) = uVm(9)+0.0001;
uSpks(9) = uSpks(9)+0.1;

fig2 = figure(); set(gcf,'color','w','Position', [50, 50, 1200, 400]);
sb5 = axes('Position',[.73,.79,.25,.16]);
left_color = [0 0 0]; right_color = [.45 .45 .45];
set(sb5,'ColorOrder',[left_color; right_color]);
yyaxis left;
xx = all_Nb(aux_ind);
semilogx(xx, uVm*1000-0.2,'*k','LineWidth',3);
ylim([-53, -51.4]);
yticks([-53, -52, -51]); ylabel('Vm mean [mV]');
yyaxis right;
semilogx(all_Nb, GcC,'color',[.75,.75,.75],'LineWidth',2); box off;
xlim([0,160000]);
xlabel('N_{ext}'); ylabel('G^{com} [%]');

sb6 = axes('Position',[.73,.56,.25,.16]);
set(sb6,'ColorOrder',[left_color; right_color]);
yyaxis left;
semilogx(xx, uSpks,'*k','LineWidth',3); hold on;
ylim([.7, 5.4]); yticks([1, 3, 5]);
ylabel('E[Spike count]');
yyaxis right;
semilogx(all_Nb, GcC,'color',[.75,.75,.75],'LineWidth',2); box off;
xlim([0,160000]); 
xlabel('N_{ext}'); ylabel('G^{com} [%]');



%Spike count distances for varying Rext -----------------------------------
load('Sim_160000_05\Spkc_Distances.mat');

sb4 = axes('Position',[.07,.08,.25,.36]);      
X = [0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13];
area(X,[cnts_s0(X(2:end)),cnts_s0(end)],'FaceColor',[.25,.25,.25],'FaceAlpha',.8); hold on;
area(X,[cnts_s3(X(2:end)),cnts_s3(end)],'FaceColor',[.9,.9,.9],'FaceAlpha',.8); 
stairs(edges,[cnts_s3,cnts_s3(end)],'color',[.6,.6,.6],'Linewidth',2);
stairs(edges,[cnts_s0,cnts_s0(end)],'color',[0,0,0],'Linewidth',2);  box off;
xlim([0,13]); ylim([0,0.4]);
xlabel('Spike count'); ylabel('P(Spike count)');
leg3 = legend({'Sparse';'Dense'}); legend boxoff;
text(4.5,.35,['1-BC = ',num2str(1-BC(8,1),2)],...
     'FontSize',13,'FontName','Arial','Fontweight','Bold');
sb4.XTick = .5:2:12.5;
sb4.XTickLabel = num2str([0:2:12]');



sb8 = axes('Position',[.45,.08,.22,.35]);
plot(1:10,1-BC(:,2),'-o','Linewidth',2,'color',[.75,.75,.75]); hold on;
plot(1:10,1-BC(:,1),'-o','Linewidth',2,'color',[.5,.5,.5]);
plot(1:10,1-BC(:,3),'-o','Linewidth',2,'color',[0,0,0]); box off;
%plot([0,11],[.1,.1],'--k'); 
%plot([0,11],[.3,.3],'--k'); box off;
xlim([0,11]);
%text(11,.05,'small'); text(11,.2,'medium'); text(11,.55,'large');
xticks([1, 5, 10]);
ylabel('1-BC'); xlabel('r_{\delta} [% of r_{bas}]');
leg4 = legend({'Smooth | Basal','Sparse | Basal','Sparse | Smooth'}); legend boxoff;


set(sb1,'Position',[.48,.15,.215,.68]); %PSTH
set(sb2,'Position',[.76,.15,.215,.68]); %Mean membrane potential
%set(sb3,'Position',[.76,.15,.215,.32]); %Std membrane potential

set(sbcond(1),'Position',[.05,.17,.18,.75]); %Conductances 1
set(sbcond(2),'Position',[.29,.17,.18,.75]); %Conductances 2
set(sbcond(3),'Position',[.53,.17,.18,.75]); %Conductances 3
set(sbcond(4),'Position',[.77,.17,.18,.75]); %Conductances 4

set(sb4,'Position',[.42,.15,.22,.76]); %Spike count distributions
set(sb5,'Position',[.07,.15,.22,.32]); %Gain and Vm
set(sb6,'Position',[.07,.62,.22,.32]); %Gain and Spk count
set(sb8,'Position',[.71,.15,.22,.76]); %Effect size

% schem = axes('Position',[.08,.60,.30,.40]);
% imshow('Figures/schematic_big.jpg');
% set(schem,'Position',[.03,.55,.40,.40]);

set(findall(fig0,'-property','FontSize'),'FontSize',12);
set(findall(fig0,'-property','FontName'),'FontName','Arial');
set(findall(fig1,'-property','FontSize'),'FontSize',12);
set(findall(fig1,'-property','FontName'),'FontName','Arial');
set(findall(fig2,'-property','FontSize'),'FontSize',12);
set(findall(fig2,'-property','FontName'),'FontName','Arial');

set(leg1,'FontSize',12,'Position',[0.485 0.85 0.08 0.06]);
set(lcond,'FontSize',12,'Position',[0.9 0.7 0.08 0.06]);
set(leg2,'FontSize',9,'Position',[0.86 0.2 0.12 0.09]);
set(leg3,'FontSize',12,'Position',[0.555 0.6 0.08 0.06]);
set(leg4,'FontSize',12,'Position',[0.71 0.81 0.13 0.09]);


% print(fig0,'fig_simulation_0','-dsvg');
% print(fig1,'fig_simulation_1','-dsvg');
% print(fig2,'fig_simulation_2','-dsvg');


