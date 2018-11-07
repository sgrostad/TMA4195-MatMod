function u = getXVelocity(x, z, kappa, m, h, d)
    
    % Velocity in x direction
    u = kappa/(m+1) * ((h(x)-d(x)).^(m+1) - (h(x)-z).^(m+1));
    
end





% function [u, v] = getVelocity(z, kappa, m, h, dhdx, d, dddx)
%     
%     % Velocity in x direction
%     u = kappa/(m+1) * ((h-d).^(m+1) - (h-z).^(m+1));
%     
%     % Velocity in z direction
%     v = kappa * ( (d-z) .* (h-d).^m .* (dhdx - dddx) + ...
%         dhdx/(m+1) * ((h-d).^(m+1) - (h-z).^(m+1)) );
%     
% end


% function [u, v] = getVelocity(z, x, kappa, m, h, dhdx, d, dddx)
%     
%     % Velocity in x direction
%     u = kappa/(m+1) * ((h(x)-d(x)).^(m+1) - (h(x)-z).^(m+1));
%     
%     % Velocity in z direction
%     v = kappa * ( (d(x)-z) .* (h(x)-d(x)).^m .* (dhdx(x) - dddx(x)) + ...
%         dhdx(x)/(m+1) * ((h(x)-d(x)).^(m+1) - (h(x)-z(x)).^(m+1)) );
%     
% end



