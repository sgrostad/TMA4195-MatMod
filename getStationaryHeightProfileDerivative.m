function dhdx = getStationaryHeightProfileDerivative(x, q, intq, d, dddx, m, kappa, J0, rho)
    % N  = Number of discretization points
    % x0 = Starting point of the model of the glacier
    % xS = Starting point snow smelts in summer
    % q0 = Accumulation rate abouve the snow line
    % a  = Slope of accumulation rate
    
    % Get adjusted height profile
    eta = ((m+2)/kappa*(intq(x) + J0/rho)).^(1/(m+2));
    
    % ... and it's derivative
    dhdx = q(x)./(kappa * eta.^(m+1)) + dddx(x);
    
end