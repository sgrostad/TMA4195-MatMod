function trajectories = getStationaryTrajectories(u, v, h, d, dddx, x0, m, kappa, rho, dt, dj, J0)
    
    % Construct uniformy spaced initial values
    niter = 10;
    z0vec = equiFlux(h, d, x0, m, kappa, rho, dj, J0, niter);
    
    
    trajectories = cell(1, length(z0vec));
    for i = 1:length(z0vec)
        [x, z] = getTrajectory(x0, z0vec(i), u, v, h, dddx, dt);
        trajectories{i} = [x; z];
    end
    

end


function [z, j] = equiFlux(h, d, x0, m, kappa, rho, dj, J0, niter) %(z0vec, u, x0, rho, j0vec, niter)

    j0vec = 0:dj:J0;
    if j0vec(end) ~= J0
        j0vec(end+1) = J0;
    end
    
    z = j0vec/J0 * (h(0) - d(0)) + d(0);
    j = zeros(size(z));
    
    for i = 1:niter
        
        j(2:end) = rho*kappa/(m+1) * ((h(0) - d(0)).^(m+1) .* (z(2:end)-d(0)) + 1/(m+1) * ((h(0) - z(2:end)).^(m+2) - (h(0) - d(0)).^(m+2)));
        j(end) = j0vec(end);
        
        z = interp1(j, z, j0vec);
        
    end
    
    z = z(2:end-1);
    x = x0 + zeros(size(z));
    
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
        
        x(end+1) = xt + u(xt,zt)*dt;
        z(end+1) = zt + v(xt,zt)*dt;
        
        if (x(end-1) == x(end)) && (z(end-1) == z(end))
            break
        end
        
    end
    
    a = 0;
    b = dt;
    for i = 1:100
        
        dt = (a+b)/2;
        
        xt = x(end-1);
        zt = z(end-1);
        
        x(end) = xt + u(xt,zt)*dt;
        z(end) = zt + v(xt,zt)*dt;
        
        if h(x(end)) >= z
            a = dt;
        else
            b = dt;
        end
    end
end