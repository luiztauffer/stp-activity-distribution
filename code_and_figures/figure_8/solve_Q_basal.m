% OUTPUT:
% PRR    - Proportion of Resources Released
% INPUT:
% tauf   - Facilitatory time constant [seconds]
% taud   - Recovery time constant [seconds]
% U      - Utilization factor
% rb     - Rate of basal synapses [spks/second]
% N      - Total number of synpases (populatoin size)
% Tb     - Burst duration [seconds]
% Re     - array of extra input rates [spks/second]

function Q = solve_Q_basal(tauf, taud, U, rb, Tb) 
    %Basal activity
    params0 = [tauf, taud, U, rb];
    y0 = init_cond(params0);
    up0 = y0(1)*(1-U)+U;
    Q = up0*y0(2)*Tb*rb;
    
    



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
