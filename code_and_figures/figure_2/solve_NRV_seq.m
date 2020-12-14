function [y, upb, Rseq] = solve_NRV_seq(tauf, taud, U, Rn, Re, N, Nb, Tb) 
    %Noisy synapses
    params0 = [tauf, taud, U, Rn];
    y0 = init_cond(params0);
    up0 = y0(1)*(1-U)+U;
    Vn = up0*y0(2)*Tb*Rn;
    
    %Bursty synapses
    tstart = 0;
    tfinal = Tb;
    Rb = Re/Nb;
    %paramsb = [tauf, taud, U, Rn+Rb];
    %[t,y] = ode45(@(t,y)myODE(t,y,paramsb),[tstart tfinal],y0);
    %y = y';
    dt = 0.0001;
    Rseq = [Rn*ones(1,1+0.5*Tb/dt),(Rb+Rn)*ones(1,Tb/dt),Rn*ones(1,0.5*Tb/dt)];
    y = evolve_cond(tauf, taud, U, Rseq, Tb, y0(1), y0(2));
    
    upb = y(1,:)*(1-U) + U;
    Vb = mean(upb.*y(2,:))*Tb*(Rn+Rb);
    
    %Combined NRV
    NRV = (N-Nb)*Vn + Nb*Vb;


    function y0 = init_cond(params)
        tauf = params(1);
        taud = params(2);
        U = params(3);
        R = params(4);

        y0(1) = U*tauf*R/(1+U*tauf*R);
        up0 = y0(1)*(1-U)+U;
        y0(2) = 1/(1+up0*taud*R);
        
        
    function y = evolve_cond(tauf, taud, U, R, T, u0, x0)

        dt = 0.0001;
        y = zeros(2,2*int32(1+T/dt)-1);
        y(1,1) = u0;
        y(2,1) = x0;
        for i = 2:2*int32(1+T/dt)-1
            y(1,i) = y(1,i-1) + dt*( -y(1,i-1)/tauf + U*(1-y(1,i-1))*R(i) );
            up = y(1,i) + U*(1-y(1,i));
            y(2,i) = y(2,i-1) + dt*( (1-y(2,i-1))/taud - up*y(2,i-1)*R(i) );
        end


    function dydt = myODE(t,y,params) %Working, but too slow
        %Y(1) = u(t)
        %Y(2) = x(t)
        
        tauf = params(1);
        taud = params(2);
        U = params(3);
        R = params(4);

        dudt = -y(1)/tauf + U*(1-y(1))*R;
        uplus = y(1)*(1-U) + U;
        dxdt = (1-y(2))/taud - uplus*y(2)*R;
        
        dydt = [dudt; dxdt];
