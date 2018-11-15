function u = getXVelocity(x, z, kappa, m, h, d)
    % Compute the velocity in x direction
    
    u = kappa/(m+1) * ((h(x)-d(x)).^(m+1) - abs(h(x)-z).^(m+1));
    
end


