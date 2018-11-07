function trajectories = getStationaryTrajectories(u, v, h, d, x0, m, kappa, rho, dt, dj)
    
    % Construct uniformy spaced initial values
    h0 = h(0);
    d0 = d(0);
    eta0 = h0 - d0;
    
    du  = dj/rho;
    umax = kappa/(m+1)*eta0^(m+1);
    uvec = du:du:umax;
    z0vec = h0 - (eta0^(m+1) - uvec*(m+1)/kappa).^(1/(m+1));
    
    trajectories = cell(1, size(uvec,1));
    for i = 1:size(uvec,1)
        
        [x, z] = getTrajectory(x0, z0vec(i), u, v, h, d, dt);
        trajectories{i} = [x; z];
        
    end

end

function [x,z] = getTrajectory(x0, z0, u, v, h, d, dt)
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
        
        x(end+1) = x(end) + u(x(end),z(end))*dt;
        z(end+1) = z(end) + v(x(end),z(end))*dt;
        
        if (x(end-1) == x(end)) && (z(end-1) == z(end))
            break
        end
        
    end

end