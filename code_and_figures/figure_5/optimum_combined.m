
Gain2 = zeros(3,length(all_noise));
Freq2 = zeros(3,length(all_noise));

i = 0;
for Rn = all_noise
    i = i+1;
    for syntype = 1:2
        if syntype == 1  %Input to target --------
            tauf = .2;
            taud = .05;
            U = .1;    
            syntitle = 's1 - facilitatory';        
        elseif syntype == 2  %Input to Inhibitory neuron --------
            tauf = .05;
            taud = .2;
            U = .7;
            syntitle = 's2 - depressive';
        end

        % Cerebellum - Parallel fibers -> Purkinje cells
        %Rn = 4.5; 
        N = 160000; 
        Tb = 0.04;

    %     % Hipoccampus - CA3 -> CA1
    %     Rn = 5.; 
    %     N = 1000; 
    %     Tb = 0.200;


        all_Nb = unique( logspace(0,4,1000) );
        all_Re = linspace(2000,8000,5);

        NRVbasal = solve_NRV(tauf, taud, U, Rn, 0, N, 1, Tb);

        nRows = length(all_Nb);
        nCols = length(all_Re);
        NRV = zeros(nRows,nCols);
        k = 1;
        for Nb = all_Nb
            j = 1;
            for Re = all_Re
                NRV(k,j) = solve_NRV(tauf, taud, U, Rn, Re, N, Nb, Tb);
                j = j+1;
            end
            k = k+1;
        end
        j = 1;
        for Re = all_Re
            NRVnoise(j) = solve_NRV(tauf, taud, U, Rn, Re, N, N, Tb)-NRVbasal;
            j = j+1;
        end

        NRV = NRV - NRVbasal;

        if syntype == 1
            NRVe_abs = NRV;        
        elseif syntype == 2  
            NRVi_abs = NRV;
        end


        %Finds Nb that maximizes NRV
        for j = 1:nCols
            [val,ind] = max( NRV(:,j) );
            nBn_max(j) = ind;

            NRV(:,j) = 100*( NRV(:,j) - NRVnoise(j)  )/NRVnoise(j) ;
        end

        if syntype == 1
            NRVe = NRV;  
            nB_maxe = nBn_max;
        elseif syntype == 2  
            NRVi = NRV;
            nB_maxi = nBn_max;
        end

    end


    %Finds Nb that maximizes combined FFI NRV
    NRVc = NRVe - kEI*NRVi;

    for j = 1:nCols
        [val,ind] = max( NRVc(:,j) );
        nBn_max(j) = ind;

        optvaln_ce(j) = NRVe(ind,j);
        optabsn_ce(j) = NRVe_abs(ind,j);
        optabsn_ee(j) = NRVe_abs(nB_maxe(j),j);
        scattern_ee(j) = NRVe_abs(end,j);

        optvaln_ci(j) = NRVi(ind,j);
        optabsn_ci(j) = NRVi_abs(ind,j);
        optabsn_ii(j) = NRVi_abs(nB_maxi(j),j);
        scattern_ii(j) = NRVi_abs(end,j);

        optvaln_cc(j) = NRVc(ind,j);
    end

    Gain2(1,i) = mean(optvaln_ce);
    Gain2(2,i) = mean(optvaln_ci);
    Gain2(3,i) = mean(optvaln_cc);
    
    %To use in LIF simulation
    optNbce = all_Nb(nB_maxe);
    optNbci = all_Nb(nB_maxi);
    optNbc = all_Nb(nBn_max);
    Freq2(1,i) = median( all_Re./optNbce );
    Freq2(2,i) = median( all_Re./optNbci );
    Freq2(3,i) = median( all_Re./optNbc ); 
    
    prr_basal = NRVbasal;
    %save('inhibition_05Hz.mat','all_Re','prr_basal','optNbc','optabs_ci','optabs_ii','scatter_ii');
end





