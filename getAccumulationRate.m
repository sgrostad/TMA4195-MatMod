function [q, xF] = getAccumulationRate(x, x0, xS, xF, q0, a)
    % Compute the accumulation rate

    % x  = Discretization points
    % x0 = Starting point of the model of the glacier
    % xS = Starting point snow smelts in summer
    % q0 = Accumulation rate abouve the snow line
    % a  = Slope of accumulation rate
    
    % Compute the toe of the glacier
    q = zeros(size(x));
    q(x>=x0 & x<xS) = q0;
    q(x>=xS & x<xF) = q0 - a*(x(x>=xS & x<xF) - xS);

end


