function v = getZVelocity(x, z, kappa, m, h, d, dhdx, dddx)
    
    % Velocity in z direction
    v = kappa * ( (d(x)-z) .* (h(x)-d(x)).^m .* (dhdx(x) - dddx(x)) + ...
        dhdx(x)/(m+1) * ((h(x)-d(x)).^(m+1) - abs(h(x)-z).^(m+1)) );
    
end
