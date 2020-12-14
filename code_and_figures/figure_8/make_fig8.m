%Generate PRR for several Gamma distributed basal rates
%Multiply synaptic PRR profile by the rate-scaled Gamma pdf
%For r_b = 0:100
clear all; close all;

% Facilitation 
tauf1 = .2;
taud1 = .05;
U1 = .1;   

% Parameters
Tb = 0.04;
all_rates_basal = linspace(0.001, 5, 1000);
rb_mean = .5;
rd = rb_mean*0.1;
[r_opt, optgain] = theoretical_optfreq(tauf1, taud1, U1, rb_mean, Tb);
OD = rd/r_opt;

% Caclulate single synapse basal/extra Q profile
i = 0;
for rb = all_rates_basal
    i = i + 1;
    Qs1_b(i) = solve_Q_basal(tauf1, taud1, U1, rb, Tb);
    Qs1_e_sparse(i) = solve_Q_extra(tauf1, taud1, U1, rb, Tb, r_opt);
    Qs1_e_dense(i) = solve_Q_extra(tauf1, taud1, U1, rb, Tb, rd);
end

%Reference values for fixed r_b=0.5
Qs1_b_r = solve_Q_basal(tauf1, taud1, U1, rb_mean, Tb);
Qs1_e_sparse_r = solve_Q_extra(tauf1, taud1, U1, rb_mean, Tb, r_opt);
Qs1_e_dense_r = solve_Q_extra(tauf1, taud1, U1, rb_mean, Tb, rd);

%Not selected units integrated expected PRR
GammaShape = linspace(1, 20, 1000); %[.5, 1, 5, 10, 20]; %shape parameter
GammaScale = rb_mean./GammaShape; %scale parameter
pdf_rates = zeros(length(GammaScale), length(all_rates_basal));
EQs1_b = zeros(1, length(GammaScale));
EQs1_e_sparse = zeros(1, length(GammaScale));
EQs1_e_dense = zeros(1, length(GammaScale));
for i = 1:length(GammaShape)
    A = GammaShape(i);
    B = GammaScale(i);
    pdf = gampdf(all_rates_basal, A, B);
    pdf_rates(i, :) = pdf/sum(pdf);   %Make sum(pdf)=1
    %rates_pdf = rates_pdf./all_rates_basal;       %Scale by rate intensity

    EQs1_b(i) = sum( pdf_rates(i,:).*Qs1_b );
    EQs1_e_sparse(i) = sum( pdf_rates(i,:).*Qs1_e_sparse );
    EQs1_e_dense(i) = sum( pdf_rates(i,:).*Qs1_e_dense );
end

% Estimate gains (array with same size as GammaShape)
EQp_b = EQs1_b;
EQp_e_dense = EQs1_e_dense;
EQp_e_sparse = EQs1_b*(1-OD) + EQs1_e_sparse*OD;
G = 100*(EQp_e_sparse-EQp_b)./(EQp_e_dense-EQp_b) - 100;

%Reference values for fixed r_b=0.5
Qp_e_sparse_r = Qs1_b_r*(1-OD) + Qs1_e_sparse_r*OD;
G_r = 100*(Qp_e_sparse_r-Qs1_b_r)./(Qs1_e_dense_r-Qs1_b_r) - 100;

% Plot figures
figure(); set(gcf,'color','w','Position', [50, 50, 1000, 410]);
sb1 = subplot(1,3,1);
plot(all_rates_basal, pdf_rates(1,:), 'k', 'LineWidth', 2); hold on; box off;
plot(all_rates_basal, pdf_rates(1000,:), 'color', [.5, .5, .5], 'LineWidth', 2);
plot([0.5, 0.5], [0, 1], '--k', 'LineWidth', .8);
xlim([0, 2]); ylim([0, 1.1*max(pdf_rates(1000,:))]);
ylabel('P(r_{bas})'); xlabel('r_{bas} [Hz]');
l1 = legend({'$$Shape = 1$$','$$Shape = 20$$', '$$\overline{r_{bas}} = 0.5 Hz$$'}); 
set(l1, 'Interpreter','latex','fontsize',12,'FontName','Arial')
legend boxoff;

sb2 = subplot(3,3,2);
plot(GammaShape, EQs1_b, 'k', 'LineWidth', 2); hold on; box off;
plt = plot([GammaShape(1), GammaShape(end)], [Qs1_b_r, Qs1_b_r], '--k', 'LineWidth', 0.8);
xlim([1, 20]);
ylabel('E[Q^s_{bas}]');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
l2 = legend(plt, {'$$fixed$$ $$r_{bas}$$'}); 
set(l2, 'Interpreter','latex','fontsize',12,'FontName','Arial'); legend boxoff;

sb3 = subplot(3,3,5);
plot(GammaShape, EQs1_e_dense, 'k', 'LineWidth', 2); hold on; box off;
plot([GammaShape(1), GammaShape(end)], [Qs1_e_dense_r, Qs1_e_dense_r], '--k', 'LineWidth', 0.8);
xlim([1, 20]);
ylabel('E[Q^s_{\delta}]');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);

sb4 = subplot(3,3,8);
plot(GammaShape, EQp_e_sparse, 'k', 'LineWidth', 2); hold on; box off;
plot([GammaShape(1), GammaShape(end)], [Qp_e_sparse_r, Qp_e_sparse_r], '--k', 'LineWidth', 0.8);
xlim([1, 20]);
xlabel('Gamma Shape');
ylabel('E[Q^s_{opt}]');


% Plot Gains estimates
sb5 = subplot(1,3,3);
flist = {'rbas_05.mat'; 'rbas_2.mat'; 'rbas_5.mat'};
colors = [0, 0.45, 0.7];
for i = 1:3
    load(char(flist(i)));
    c = colors(i);
    plts(i) = plot(GammaShape, G, 'color', [c, c, c], 'LineWidth', 2); hold on;
    plts2(i) = plot([GammaShape(1), GammaShape(end)], [G_r, G_r], '--', 'color', [c, c, c], 'LineWidth', 0.8);
    xlim([1, 20]); ylim([0, 70]);
    xlabel('Gamma Shape'); ylabel('G_{opt} [%]');
end
box off;
l5 = legend([plts, plts2(1)], {'$$\overline{r_{bas}} = 0.5 Hz$$','$$\overline{r_{bas}} = 2 Hz$$','$$\overline{r_{bas}} = 5 Hz$$', '$$fixed$$ $$r_{bas}$$'}); 
set(l5, 'Interpreter','latex','fontsize',12,'FontName','Arial','FontAngle','normal'); legend boxoff;

fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',12);
set(findall(fig,'-property','FontName'),'FontName','Arial');

sb1.Position = [.10, .14, .213, .79];
sb2.Position = [0.41, 0.72, 0.213, 0.21];
sb3.Position = [0.41, 0.43, 0.213, 0.21];
sb4.Position = [0.41, 0.14, 0.213, 0.21];
sb5.Position = [0.72, 0.14, 0.213, 0.79];

%save('rbas_5.mat', 'G', 'G_r');
