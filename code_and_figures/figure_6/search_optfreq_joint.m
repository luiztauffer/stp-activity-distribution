function [Freq, Gain] = search_optfreq_joint(tauf1, taud1, U1, tauf2, taud2, U2,...
                                                   Rn, Tb)


for syntype = 1:2
    if syntype == 1  %Input to target --------
        tauf = tauf1;
        taud = taud1;
        U = U1;           
    elseif syntype == 2  %Input to Inhibitory neuron --------
        tauf = tauf2;
        taud = taud2;
        U = U2;
    end

    % Cerebellum - Parallel fibers -> Purkinje cells
    %N = 160000; 
    N = 10000; 

    all_Nb = unique( logspace(0,5,1200) );
    all_Re = linspace(4000,8000,3);

    NRVbasal = solve_NRV(tauf, taud, U, Rn, 0, N, 1, Tb);

    nRows = length(all_Nb);
    nCols = length(all_Re);
    NRV = zeros(nRows,nCols);
    k = 1;
    for Nb = all_Nb
        j = 1;
        for Re = all_Re
            NRV(k,j) = solve_NRV(tauf, taud, U, Rn, Re, N, Nb, Tb) - NRVbasal;
            j = j+1;
        end
        k = k+1;
    end
    j = 1;
    for Re = all_Re
        NRVnoise(j) = solve_NRV(tauf, taud, U, Rn, Re, N, N, Tb) - NRVbasal;
        j = j+1;
    end

    if syntype == 1
        NRVe_abs = NRV;        
    elseif syntype == 2  
        NRVi_abs = NRV;
    end

    %Finds Nb that maximizes NRV
    for j = 1:nCols
        [val,ind] = max( NRV(:,j) );
        nB_max(j) = ind;

        NRV(:,j) = 100*( NRV(:,j) - NRVnoise(j)  )/NRVnoise(j) ;
    end

    if syntype == 1
        NRVe = NRV;  
        nB_maxe = nB_max;
    elseif syntype == 2  
        NRVi = NRV;
        nB_maxi = nB_max;
    end

end


%Finds Nb that maximizes combined FFI NRV
NRVc = NRVe - NRVi;

for j = 1:nCols
    [val,ind] = max( NRVc(:,j) );
    nB_max(j) = ind;

    optval_ce(j) = NRVe(ind,j);
    optabs_ce(j) = NRVe_abs(ind,j);
    optabs_ee(j) = NRVe_abs(nB_maxe(j),j);
    distri_ee(j) = NRVe_abs(end,j);

    optval_ci(j) = NRVi(ind,j);
    optabs_ci(j) = NRVi_abs(ind,j);
    optabs_ii(j) = NRVi_abs(nB_maxi(j),j);
    distri_ii(j) = NRVi_abs(end,j);

    optval_cc(j) = NRVc(ind,j);
end

%Gain = zeros(1,3);
%Freq = zeros(1,3);

%Gain(1) = mean(optval_ce);
%Gain(2) = mean(optval_ci);
Gain = mean(optval_cc);

%To use in LIF simulation
optNbce = round(all_Nb(nB_maxe));
optNbci = round(all_Nb(nB_maxi));
optNbc = round(all_Nb(nB_max));
%Freq(1) = median( all_Re./optNbce );
%Freq(2) = median( all_Re./optNbci );
Freq = median( all_Re./optNbc ); 

prr_basal = NRVbasal;


