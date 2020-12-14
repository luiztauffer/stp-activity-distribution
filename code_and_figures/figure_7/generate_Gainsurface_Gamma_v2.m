%Generate PRR and Gains for several Gamma distributions
%Multiply synaptic PRR profile by the rate-scaled Gamma pdf
%For rate = delta:8000
%E[PRRp] = sum( Gamma(A,B,rate)*PRRs(rate)/rate )
clear all; close all;

rb = .5;  
Tb = 0.04;
N = 160000;
Re = 8000;

%Facilitation --------
tauf1 = .2;
taud1 = .05;
U1 = .1;   
%Depression --------
tauf2 = .05;
taud2 = .2;
U2 = .7;

rate_delta = Re/N;
all_rates = linspace(rate_delta,Re,50000);

%Caclulate single synapse basal PRR
BasalPRR1 = solve_PRR_modelfree(tauf1, taud1, U1, rb, 1, Tb, 0);
BasalPRR2 = solve_PRR_modelfree(tauf2, taud2, U2, rb, 1, Tb, 0);

%Caclulate single synapse PRR profile (per extra rate)
i = 0;
for rate = all_rates
    i = i + 1;
    PRRs1(i) = solve_PRR_modelfree(tauf1, taud1, U1, rb, 1, Tb, rate) - BasalPRR1;
    PRRs2(i) = solve_PRR_modelfree(tauf2, taud2, U2, rb, 1, Tb, rate) - BasalPRR2;
end

%Caclulate synapse PRR for smooth distribution (rate_delta increase)
SmoothPRRp1 = 1*PRRs1(1)/rate_delta;
SmoothPRRp2 = 1*PRRs2(1)/rate_delta;

GammaShape = linspace(1, 20, 153); %[.5, 1, 5, 10, 20]; %shape parameter
GammaScale = logspace(-1.5, 3, 600); %[.1, .5, 1, 2, 5, 10, 100, 200, 300, 500, 1000]; %scale parameter
PRRp1 = zeros(length(GammaShape),length(GammaScale));
PRRp2 = zeros(length(GammaShape),length(GammaScale));
Next = zeros(length(GammaShape),length(GammaScale));
i = 0;
for A = GammaShape
    i = i + 1;
    j = 0;
    for B = GammaScale
        j = j + 1;
        rates_pdf = gampdf(all_rates,A,B);
        rates_pdf = rates_pdf/sum(rates_pdf);   %Make sum(pdf)=1
        rates_pdf = rates_pdf./all_rates;       %Scale by rate intensity

        PRRp1(i,j) = sum( rates_pdf.*PRRs1 );
        PRRp2(i,j) = sum( rates_pdf.*PRRs2 );
    end
end

Gain1 = 100*(PRRp1 - SmoothPRRp1)/SmoothPRRp1;
Gain2 = 100*(PRRp2 - SmoothPRRp2)/SmoothPRRp2;

save('Expected_Gain.mat','PRRp1','PRRp2','Gain1','Gain2','GammaShape','GammaScale');

figure();
subplot(1,2,1);
imagesc(Gain1); set(gca,'XDir','reverse');
colormap(brewermap([],'*RdBu')); c1 = colorbar();
subplot(1,2,2);
imagesc(Gain2); set(gca,'XDir','reverse');
colormap(brewermap([],'*RdBu')); c2 = colorbar();
