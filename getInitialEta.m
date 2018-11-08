function [etaInit, xF] = getInitialEta(x0, xS, q0, a, J0, rho, m, kappa)
    xF = getStationaryToe(x0, xS, q0, a, J0, rho);

    % Accumulation rate
%     q    = @(x) getAccumulationRate(x, x0, xS, xF, q0, a);
    intq = @(x) getCumulativeAccumulationRate(x, x0, xS, xF, q0, a, J0, rho);


    % Bedrock profile
    d    = @(x) 0.02*(1+sin(pi*x/xS + 0));
%     dddx = @(x) 0.02*pi/xS * cos(pi*x/xS + 0);

    % Glacier height profile
    h    = @(x) getStationaryHeightProfile(x, intq, d, m, kappa, J0, rho);
%     dhdx = @(x) getStationaryHeightProfileDerivative(x, q, intq, d, dddx, m, kappa, J0, rho);
    
    etaInit =@(x) h(x) - d(x);
end

