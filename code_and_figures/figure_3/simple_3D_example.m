% Simple 3D example
clear all; close all;

Tb = 0.04;
R0 = 0; 
Re = 100;

for syntype = 1:3
    if syntype == 1
        %Facilitation --------
        tauf = .2;
        taud = .05;
        U = .1; 
    elseif syntype == 2
        %Depression --------
        tauf = .05;
        taud = .2;
        U = .7;
    elseif syntype == 3
        %Depression/Facilitation ----------
        tauf = 0.1;
        taud = 0.1;
        U = 0.4;   
    end
    
    i = 0;
    for rate_dist = 0:.05:.5
        i = i + 1;
        re = rate_dist*Re;
        [~, Vs] = solve_NRV(tauf, taud, U, R0, re, Tb);
        PRR(syntype,i) = Vs;
    end
end

fig1 = figure(); set(gcf,'color','w','Position', [50, 50, 600, 400]);
plot(PRR(1,:)); hold on;
plot(PRR(2,:));
plot(PRR(3,:)); box off;