function trajectories = getStationaryTrajectories(u, v, h, d, dddx, x0, m, kappa, rho, dt, dJ)
    
    % Construct uniformy spaced initial values
    h0 = h(x0);
    d0 = d(x0);
    eta0 = h0 - d0;
    
    du  = dJ/rho;
    umax = kappa/(m+1)*eta0^(m+1);
    uvec = du:du:umax;
    z0vec = h0 - (eta0^(m+1) - uvec*(m+1)/kappa).^(1/(m+1));
    
    dz = 0.05;
    z0vec = d(0)+dz:dz:h(0);
    
    trajectories = cell(1, length(z0vec));
    for i = 1:length(z0vec)
        [x, z] = getTrajectory(x0, z0vec(i), u, v, h, dddx, dt);
        trajectories{i} = [x; z];
        i
    end
    

end

function [x,z] = getTrajectory(x0, z0, u, v, h, dddx, dt)
    % Solves the pde given by
    %   x' = u(x,z)
    %   z' = v(x,z)
    % with initial condition
    %   x(0) = x0
    %   z(0) = z0
    % until z > h
    
    x = x0;
    z = z0;
    
    while h(x(end)) >= z
        
        xt = x(end);
        zt = z(end);
        
        x(end+1) = xt + u(x(end),z(end))*dt;
        z(end+1) = zt + (v(x(end),z(end)) + 0*dddx(x(end-1)))*dt;
        
        if (x(end-1) == x(end)) && (z(end-1) == z(end))
            break
        end
        
    end

end