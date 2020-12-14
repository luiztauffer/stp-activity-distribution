%Calculates the effect size for the Mann–Whitney U test with the
%Common Language effect size.
% McGraw, K. O., & Wong, S. P. (1992). A common language effect size statistic. Psychological Bulletin, 111(2), 361-365.
%"how often a score sampled from one group will be greater than a score 
%sampled from another group"
clear all; close all;

for k = 1:10
    load(['Sim_160000_05/New_Neuron_sim_1_perc_',num2str(k),'.mat']);
    spkt1 = spkt;
    VmP1 = VmP;
    load(['Sim_160000_05/New_Neuron_sim_2_perc_',num2str(k),'.mat']);
    spkt2 = spkt;
    VmP2 = VmP;
    load(['Sim_160000_05/New_Neuron_sim_3_perc_',num2str(k),'.mat']);
    spkt3 = spkt;
    VmP3 = VmP;
    load(['Sim_160000_05/New_Neuron_sim_0_perc_',num2str(k),'.mat']);
    spkt0 = spkt;
    VmP0 = VmP;

    nTrials = 3010;
    nBins = 120;
    dt = 0.1;

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
    SpkCnt_b(1:10) = [];
    SpkCnt_s0(1:10) = [];
    SpkCnt_s1(1:10) = [];
    SpkCnt_s2(1:10) = [];
    SpkCnt_s3(1:10) = [];
    nTrials = nTrials - 10;

    edges = 0:13;
    [cnts_b,~] = histcounts(SpkCnt_b,edges,'Normalization','probability');
    [cnts_s0,~] = histcounts(SpkCnt_s0,edges,'Normalization','probability');
    [cnts_s3,~] = histcounts(SpkCnt_s3,edges,'Normalization','probability');
    
    s0_b_diff = [];
    s3_b_diff = [];
    s0_s3_diff = [];
%     for i = 1:nTrials
%         bas = SpkCnt_b(i);
%         s3 = SpkCnt_s3(i);
%         s0_b_diff(i,:) = SpkCnt_s0 - bas;
%         s0_s3_diff(i,:) = SpkCnt_s0 - s3;
%         s3_b_diff(i,:) = SpkCnt_s3 - bas;
%     end
%     s0_b_diff = reshape(s0_b_diff,1,[]);
%     s0_s3_diff = reshape(s0_s3_diff,1,[]);
%     s3_b_diff = reshape(s3_b_diff,1,[]);
% 
%     CL(k,1) = sum(s0_b_diff>0)/(nTrials^2);
%     CL(k,2) = sum(s3_b_diff>0)/(nTrials^2);
%     CL(k,3) = sum(s0_s3_diff>0)/(nTrials^2);
    
    
    % Effect size for the Wilcoxon rank sum test
    [mwp(k,1),~,st] = ranksum(SpkCnt_b,SpkCnt_s0);
    MWES(k,1) = abs( st.zval/sqrt(2*nTrials) );
    [mwp(k,2),~,st] = ranksum(SpkCnt_b,SpkCnt_s3);
    MWES(k,2) = abs( st.zval/sqrt(2*nTrials) );
    [mwp(k,3),~,st] = ranksum(SpkCnt_s0,SpkCnt_s3);
    MWES(k,3) = abs( st.zval/sqrt(2*nTrials) );
    
    %Bhattacharya Coefficient
    BC(k,1) = sum(sqrt(cnts_b.*cnts_s0));
    BC(k,2) = sum(sqrt(cnts_b.*cnts_s3));
    BC(k,3) = sum(sqrt(cnts_s0.*cnts_s3));
    
    %Bhattacharya Distance
    BD(k,1) = -log( BC(k,1) );
    BD(k,2) = -log( BC(k,2) );
    BD(k,3) = -log( BC(k,3) );
    
    %Hellinger distance
    HD(k,1) = sqrt( sum((sqrt(cnts_b)-sqrt(cnts_s0)).^2) )/sqrt(2);
    HD(k,2) = sqrt( sum((sqrt(cnts_b)-sqrt(cnts_s3)).^2) )/sqrt(2);
    HD(k,3) = sqrt( sum((sqrt(cnts_s0)-sqrt(cnts_s3)).^2) )/sqrt(2);
    
end

save('Spkc_Distances.mat','HD','BD','BC','MWES');

subplot(2,2,1);
plot(1-BC(:,1)); hold on
plot(1-BC(:,2));
plot(1-BC(:,3));
title('1-Bhattacharya Coefficient');

subplot(2,2,2);
plot(MWES(:,1)); hold on
plot(MWES(:,2));
plot(MWES(:,3));
title('Mann-Whitney ES');

subplot(2,2,3);
plot(BD(:,1)); hold on
plot(BD(:,2));
plot(BD(:,3));
title('Bhattacharya Distance');

subplot(2,2,4);
plot(HD(:,1)); hold on
plot(HD(:,2));
plot(HD(:,3));
title('Hellinger Distance');
