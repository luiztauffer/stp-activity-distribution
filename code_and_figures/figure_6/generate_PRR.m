%
% Generate PRR values to be used on LIF simulation
%
%--------------------------------------------------------------------------
clear all; close all;

Tb = 0.04;

tauf1 = .2;
taud1 = .05;
U1 = 0.1;
tauf2 = .05;
taud2 = .2;
U2 = .7;
    
Rn = .5; %[.5, 1, 5];

for perc = 1:10
    Re = 80000*perc/100;

    %Network size changes to keep total Rbas = 80kHz
    N = 160000; %80000/Rn

    %Search optimum value at the G surface
    [Freq, Gain] = search_optfreq_joint(tauf1, taud1, U1, tauf2, taud2, U2, Rn, Tb);
    Nopt = round(Re/Freq);

    %absolut value for PRR basal (ongoing activity)
    PRRE(1) = solve_NRV(tauf1, taud1, U1, Rn, 0, N, 1, Tb);
    PRRI(1) = solve_NRV(tauf2, taud2, U2, Rn, 0, N, 1, Tb);

    %absolut value for PRR at optimal spatial code
    PRRE(2) = solve_NRV(tauf1, taud1, U1, Rn, Re, N, Nopt, Tb);
    PRRI(2) = solve_NRV(tauf2, taud2, U2, Rn, Re, N, Nopt, Tb);

    %absolut value for PRR at distributed spatial code (noise fluctuation)
    PRRE(3) = solve_NRV(tauf1, taud1, U1, Rn, Re, N, N, Tb);
    PRRI(3) = solve_NRV(tauf2, taud2, U2, Rn, Re, N, N, Tb);


    %Gain over noise at Excitatory branch
    GcE = 100*((PRRE-PRRE(1))/(PRRE(3)-PRRE(1)) -1);
    %Gain over noise at Inhibitory branch
    GcI = 100*((PRRI-PRRI(1))/(PRRI(3)-PRRI(1)) -1);
    %Combined Gain over noise
    GcC = GcE - GcI;

    save(['PRR_sim_',num2str(perc),'perc.mat'],'N','Rn','Re','Nopt','PRRE','PRRI');
end


%save('PRRsim_1perc.mat','N','Rn','Re','Nopt','PRRE','PRRI');
%save('PRRsim_2perc.mat','N','Rn','Re','Nopt','PRRE','PRRI');
%save('PRRsim_3perc.mat','N','Rn','Re','Nopt','PRRE','PRRI');
%save('PRRsim_4perc.mat','N','Rn','Re','Nopt','PRRE','PRRI');
%save('PRRsim_5perc.mat','N','Rn','Re','Nopt','PRRE','PRRI');
%save('PRRsim_6perc.mat','N','Rn','Re','Nopt','PRRE','PRRI');
%save('PRRsim_7perc.mat','N','Rn','Re','Nopt','PRRE','PRRI');
%save('PRRsim_10perc.mat','N','Rn','Re','Nopt','PRRE','PRRI');
%save('PRRsim_20perc.mat','N','Rn','Re','Nopt','PRRE','PRRI');