function eta = finiteVolume(etaInit, intX, dx, L, dt, timeMax, xF ,x0, xS, q0, a, J0, m, kappa,rho)
    f = @(eta) kappa/(m+2) * eta.^(m+2);
    eta = zeros(length(intX), floor(timeMax/dt+1));
    if length(etaInit) == 1 %etaInit is a function handle
        eta(:,1) = etaInit(intX);
    else
        eta(:,1) = etaInit;
    end
    q    = @(x,xF) getAccumulationRate(x, x0, xS, xF, q0, a);
    for n = 1:timeMax/dt
        qHat = q(intX,xF);
        eta(:,n+1) = eta(:,n) - dt/dx * (f(eta(1:end,n)) - [J0/rho; f(eta(1:end-1,n))]) ...
           + dt*qHat;
        eta(eta(:,n+1)<0,n+1) = 0;
        eta(isnan(eta(:,n+1)),n+1) = 0;
        xF = (find(eta(:,n+1)<=1e-15,1) - 1)/length(intX)*L; 
        if sum(size(xF)) ~= 2 % No zero elements in eta
            warning('Eta array is full')
            xF = L;
        end
    end
end