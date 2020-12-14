function [optfreq, optgain] = theoretical_optfreq_joint(tauf1, taud1, U1, tauf2, taud2, U2,... 
                                             kEI, Rn, Tb) 
    %Noisy synapses
    params0 = [tauf1, taud1, U1, Rn];
    y01 = init_cond(params0);
    up01 = y01(1)*(1-U1)+U1;
    Vn1 = up01*y01(2)*Tb*Rn;
    
    params0 = [tauf2, taud2, U2, Rn];
    y02 = init_cond(params0);
    up02 = y02(1)*(1-U2)+U2;
    Vn2 = up02*y02(2)*Tb*Rn;
    
    Rext = Rn:1:800;
    NRV1 = zeros(1,length(Rext));
    NRV2 = zeros(1,length(Rext));
    i = 0;
    for Rb = Rext
        i = i + 1;
        paramsb = [tauf1, taud1, U1, Rn+Rb, Tb, y01(1), y01(2)];
        y1 = evolve_cond(paramsb);
        upb1 = y1(1,:)*(1-U1) + U1;
        NRV1(i) = mean(upb1.*y1(2,:))*Tb*(Rn+Rb);
        
        paramsb = [tauf2, taud2, U2, Rn+Rb, Tb, y02(1), y02(2)];
        y2 = evolve_cond(paramsb);
        upb2 = y2(1,:)*(1-U2) + U2;
        NRV2(i) = mean(upb2.*y2(2,:))*Tb*(Rn+Rb);
    end
    
    Re_step = Rext(2)-Rext(1);
    aux0 = Vn1 - NRV1(2:end) + Rext(2:end).*diff(NRV1)/Re_step;
    aux1 = Vn2 - NRV2(2:end) + Rext(2:end).*diff(NRV2)/Re_step;
    aux2 = kEI*(Vn1+NRV1(1))/(Vn2+NRV2(1));
    aux3 = aux2*aux1-aux0;

    zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
    zx = zci(aux3);
    if isempty(zx)
        ind = 1;
    else ind = zx(1);
    end

    optfreq = Rext(ind);
    
    S1 = Vn1/Rn;
    S2 = Vn2/Rn;
    G1 = 100*( (NRV1(ind)-Vn1)/(optfreq*S1) - 1 ); 
    G2 = 100*( (NRV2(ind)-Vn2)/(optfreq*S2) - 1 );
    
    optgain = G1 - kEI*G2;

    

    function y0 = init_cond(params)
        tauf = params(1);
        taud = params(2);
        U = params(3);
        R = params(4);

        y0(1) = U*tauf*R/(1+U*tauf*R);
        up0 = y0(1)*(1-U)+U;
        y0(2) = 1/(1+up0*taud*R);
        
        
    function y = evolve_cond(params)
        tauf = params(1);
        taud = params(2);
        U = params(3);
        R = params(4);  
        T = params(5);
        u0 = params(6);
        x0 = params(7);

        dt = 0.0001;
        y = zeros(2,int32(1+T/dt));
        y(1,1) = u0;
        y(2,1) = x0;
        for i = 2:int32(1+T/dt)
            y(1,i) = y(1,i-1) + dt*( -y(1,i-1)/tauf + U*(1-y(1,i-1))*R );
            up = y(1,i) + U*(1-y(1,i));
            y(2,i) = y(2,i-1) + dt*( (1-y(2,i-1))/taud - up*y(2,i-1)*R );
        end


