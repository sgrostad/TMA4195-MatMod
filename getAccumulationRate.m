function [q, intq, xF] = getAccumulationRate(x, x0, xS, q0, a, J0, rho)
    
    [q, intq, xF] = piecewiseLinear(x, x0, xS, q0, a, J0, rho);
    
end


function [q, intq, xF] = piecewiseLinear(x, x0, xS, q0, a, J0, rho)
    
    % x  = Discretization points
    % x0 = Starting point of the model of the glacier
    % xS = Starting point snow smelts in summer
    % q0 = Accumulation rate abouve the snow line
    % a  = Slope of accumulation rate
    
    % Compute the toe of the glacier
    xF = 1/a * (q0 + a*xS + sqrt(...
         q0^2 + 2*a*(J0/rho + q0*(xS-x0) + a*xS^2 - a*xS)));
    
    q = zeros(size(x));
    q(x>=x0 & x<xS) = q0;
    q(x>=xS & x<xF) = q0 - a*(x(x>=xS & x<xF) - xS);
    
    intq = zeros(size(x)) - J0/rho;
    intq(x>=x0 & x<xS) = q0*(x(x>=x0 & x<xS)-x0);
    intq(x>=xS & x<xF) = q0*(xS-x0) + q0*(x(x>=xS & x<xF)-xS) + a*xS*(x(x>=xS & x<xF)-xS) - a/2*(x(x>=xS & x<xF).^2 - xS.^2);

end

