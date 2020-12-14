function [optfreq, optgain] = theoretical_optfreq(tauf, taud, U, Rn, Tb) 
    %Noisy synapses
    params0 = [tauf, taud, U, Rn];
    y0 = init_cond(params0);
    up0 = y0(1)*(1-U)+U;
    Vn = up0*y0(2)*Tb*Rn;
    
    Rext = .1:.1:250;
    NRV = zeros(1,length(Rext));
    i = 0;
    for Rb = Rext
        i = i + 1;
        paramsb = [tauf, taud, U, Rn+Rb, Tb, y0(1), y0(2)];
        y = evolve_cond(paramsb);

        upb = y(1,:)*(1-U) + U;
        NRV(i) = mean(upb.*y(2,:))*Tb*(Rn+Rb);
    end
    
    Re_step = Rext(2)-Rext(1);
    aux = Re_step*( NRV(2:end)-Vn )./diff(NRV);
    %[~,ind] = min( abs( Rext(2:end) - aux ) );
    %zero-crossing point of the optimum freq equation
    zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
    zx = zci(Rext(2:end) - aux);
    if isempty(zx)
        ind = 1;
    else ind = zx(1);
    end

    optfreq = Rext(ind);
    
    Sb = Vn/Rn;
    appDnrv = NRV(Rext==optfreq) - Vn ;
    appgain = 100*( appDnrv/(optfreq*Sb) -1 );
    
    %Rd = Rn*linspace(0.01,1,2);
    Rd = Rn/10;
    NRVd = zeros(1,length(Rd));
    i = 0;
    for rd = Rd
        i = i + 1;
        paramsb = [tauf, taud, U, Rn+rd, Tb, y0(1), y0(2)];
        y = evolve_cond(paramsb);

        upb = y(1,:)*(1-U) + U;
        NRVd(i) = mean(upb.*y(2,:))*Tb*(Rn+rd);
    end
    
    Dnrv = (NRV(Rext==optfreq) - Vn)./(NRVd - Vn);
    optgain = 100*( Rd.*Dnrv/optfreq -1 );
    
    


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


