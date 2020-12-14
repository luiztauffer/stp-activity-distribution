%Calculates the PRR surface and estimates the optimum number of recruited
%fibers that maximizes the PRR of the presynaptic population given the
%Tb and Rext
%clear all; %close all;


%Depression / Facilitation --------
% tauf = .2;
% taud = .005;
% U = .2;

%Facilitation --------
tauf = .2;
taud = .05;
U = .1;

Rn = 0.5; 
N = 160000; 
Tb = 0.04;

all_Nb = unique( logspace(0,5,1000) );
all_Re = linspace(1600,8000,60);

NRVbasal = solve_NRV(tauf, taud, U, Rn, 0, N, 1, Tb);

nRows = length(all_Nb);
nCols = length(all_Re);
NRV = zeros(nRows,nCols);
j = 1;
for Re = all_Re
    i = 1;
    for Nb = all_Nb
        NRV(i,j) = solve_NRV(tauf, taud, U, Rn, Re, N, Nb, Tb);
        i = i+1;
    end
    j = j+1;
end

%NRV = 100*(NRV - NRVbasal)/NRVbasal;
NRVext = NRV - NRVbasal;
NRV = NRV - NRVbasal;

%Finds Nb that maximizes NRV
nB_max = zeros(1,nCols);
PRR_maxval = zeros(1,nCols);
for j = 1:nCols
    [val,ind] = max( NRV(:,j) );
    nB_max(j) = ind;
    PRR_maxval(j) = val;
    G_maxval(j) = 100*( val - NRV(end,j) )/NRV(end,j);

    Gain(:,j) = 100*( NRV(:,j) - NRV(end,j) )/NRV(end,j);
end


subplot('Position',[.38 .1 .23 .36]);
imagesc(Gain);
set(gca,'YDir','normal'); hold on;
plot(nB_max,'k','LineWidth',2);
[~, indx] = min(abs(all_Re-6400));
plot([indx, indx],[0 1000],'color',[.5, .5, .5],'LineWidth',2);
colormap('jet'); c = colorbar();
title(c,'G [%]');
h = gca();
h.XTick = linspace(1,nCols,5);
h.YTick = linspace(1,nRows,5);
h.XTickLabel = strread(num2str([2 4 6 8 10]),'%s');
h.YTickLabel = cellstr(num2str(round(log10([1,10,100,1000,10000]')), '10^%d'));
xlabel('r_{\delta} [% of r_{bas}]');
ylabel('N_{ext}');
set(c,'Position',[.615 .1 .02 .36]);

% % Create colormap blue for negative, red for positive and white for 0
% redColorMap = [zeros(1, 132), linspace(0, 1, 124)];
% greenColorMap = [linspace(0, 1, 124), linspace(1, 0, 132)];
% blueColorMap = [linspace(1, 0, 124), zeros(1, 132)];
% colorMap = [redColorMap; greenColorMap; blueColorMap]';
% % Apply the colormap.
% colormap(colorMap);


% subplot('Position',[.09 .1 .2 .36]);
% semilogx(all_Nb,Gain(:,indx),'k','LineWidth',2); hold on;
% semilogx(all_Nb,zeros(1,length(all_Nb)),'--k'); box off;
% [val, indy] = max( Gain(:,indx) );
% hmax = plot(all_Nb(indy), val, 'ok','LineWidth',2);
% ylabel('G [%]'); xlabel('N_{ext}');
% l2 = legend(hmax,['N_{opt} = ',num2str(round(all_Nb(indy)))]); legend boxoff;
% set(l2,'Position',[0.137 0.366 0.12 0.05]);
% xlim([0,100000]); xticks([10, 1000, 100000]); %xticks([1, 10, 100, 1000, 10000]);


subplot('Position',[.75 .33 .2 .14]);
plot(all_Re, all_Nb(nB_max),'k','LineWidth',2); box off;
ylabel('N_{opt}');
xlabel('r_{\delta} [% of r_{bas}]');
iso_freq = mean( all_Re./all_Nb(nB_max) );
xlim([1600,8000]); xticks([1600, 3200, 4800, 6400, 8000]);
h = gca(); h.XTickLabel = strread(num2str([2 4 6 8 10]),'%s');
l3 = legend(['r_{opt} = ',num2str(round(iso_freq)),' Hz']); legend boxoff;
set(l3,'Position',[0.75 0.42 0.15 0.05]);






% Extra noise distance -------------------------------------------------
extradistance = 0;
if extradistance == 1
    N = 10000;
    NRVbasal = solve_NRV(tauf, taud, U, Rn, 0, N, 1, Tb);

    all_Noise = linspace(2000,40000,381);
    NRVextnoise = zeros(1,length(all_Noise));
    j = 1;
    for Re = all_Noise
        NRVextnoise(j) = solve_NRV(tauf, taud, U, Rn, Re, N, N, Tb);
        j = j+1;
    end
    NRVextnoise = NRVextnoise - NRVbasal;
    j = 1;
    for prrval = PRR_maxval
        [val,ind] = min( abs(NRVextnoise - prrval) );
        Extra_noise(j) = all_Noise(ind);
        j = j+1;
    end

    figure();
    subplot(2,1,1);
    plot(NRVextnoise); hold on;
    plot(PRR_maxval);
    subplot(2,1,2);
    plot(all_Re,Extra_noise);
end

    
%Log-normal calc ------------------------------------------------------
lognorm = 0;
if lognorm== 1
    x = all_Nb;
    y = Gain(:,indx);
    %x(y<0) = [];
    y(y<0) = 0;
    y = y/sum(y);
    %y = circshift(y,13);
    %x = x-x(1);

    data = randsample(x,50*length(x),true,y);
    xx = log(data);
    phat = binofit(xx,1000);
    % [par,parmci] = lognfit(data);
    % Y = lognpdf(x,par(1),par(2));
    % Y = Y/sum(Y);
    par = wblfit(data);
    Y = wblpdf(x,par(1),par(2));
    Y = Y/sum(Y);
    % par = mle(xx,'pdf',@(xx,u,var) pdf('norm',xx,u,var),'start',[2,2] );
    % Y = normpdf(log(x),par(1),par(2));
    % Y = Y/sum(Y);
    %par = mle(xx,'distribution','stable');
    %pd_levy = makedist('Stable','alpha',par(1),'beta',par(2),'gam',par(3),'delta',par(4));
    %Y = pdf(pd_levy,x);
    %par = mle(xx,'pdf',@(xx,gam) pdf('Stable',xx,.5,1,gam,15),'start',2 );
    %pd_levy = makedist('Stable','alpha',0.5,'beta',1,'gam',par,'delta',-13);
    %Y = pdf(pd_levy,x);

    figure();
    subplot(2,1,1);
    plot(x,y); hold on;
    plot(x,Y);
    subplot(2,1,2);
    plot(log(x),y); hold on;
    plot(log(x),Y);
end






