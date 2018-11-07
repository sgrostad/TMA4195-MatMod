function [u, v] = getVelocity(z, kappa, m, h, dhdx, d, dddx)
    
    % Velocity in x direction
    u = kappa/(m+1) * ((h-d).^(m+1) - (h-z).^(m+1));
    
    % Velocity in z direction
    v = kappa * ( (d-z) .* (h-d).^m .* (dhdx - dddx) + ...
        dhdx/(m+1) * ((h-d).^(m+1) - (h-z).^(m+1)) );
    
end

