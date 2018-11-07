function [h, dhdx] = getStationaryHeightProfile(x, x0, xS, q0, a, m, kappa, J0, rho, d, dddx)
    % N  = Number of discretization points
    % x0 = Starting point of the model of the glacier
    % xS = Starting point snow smelts in summer
    % q0 = Accumulation rate abouve the snow line
    % a  = Slope of accumulation rate
    
    % Compute the toe of the glacier
    [q, intq] = getAccumulationRate(x, x0, xS, q0, a, J0, rho);

    eta = ((m+2)/kappa*(intq + J0/rho)).^(1/(m+2));
    
    h = eta + d;
    dhdx = q./(kappa * eta.^(m+1)) + dddx;
    
end
