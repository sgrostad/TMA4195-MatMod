function [trajectories, times] = getStationaryTrajectories(u, v, h, d, dddx, intq, x0, m, kappa, rho, dt, dj, J0, maxxqpos)
    
    % Construct uniformy spaced initial values
    niter = 10;
    [z0vec, x0vec] = equiFlux(h, d, intq, x0, m, kappa, rho, dj, J0, maxxqpos, niter);
    trajectories = cell(1, length(z0vec));
    times = zeros(1, length(z0vec));
    for i = 1:length(z0vec)
        [x, z, t] = getTrajectory(x0vec(i), z0vec(i), u, v, h, dddx, dt);
        trajectories{i} = [x; z];
        times(i) = t;
        plot(x, z, 'k');
    end
    
end



function [z, x] = equiFlux(h, d, intq, x0, m, kappa, rho, dj, J0, maxxqpos, niter) %(z0vec, u, x0, rho, j0vec, niter)
    
    % Estimate the z-values for x = x0 that gives constant flux between
    % the x's by fixed point iteration.
    
    j0vec = 0:dj:J0;
    if j0vec(end) ~= J0
        j0vec(end+1) = J0;
    end
    
    zz = j0vec/J0 * (h(0) - d(0)) + d(0);
    j = zeros(size(zz));
    
    for i = 1:niter
        j(2:end) = rho*kappa/(m+1) * ((h(0) - d(0)).^(m+1) .* (zz(2:end)-d(0)) + 1/(m+1) * ((h(0) - zz(2:end)).^(m+2) - (h(0) - d(0)).^(m+2)));
        j(end) = j0vec(end);
        zz = interp1(j, zz, j0vec);
    end
    
    z = zz(2:end-1);
    x = x0 + zeros(size(z));
    
    % Estimate the x's for z = h giving constant flux between the points.
    
    qtot = intq(maxxqpos);
    j0vec = (j0vec(end-1)+dj-J0):dj:(rho*qtot);
    if isempty(j0vec)
        return
    end
    if j0vec(end) ~= rho*qtot
        j0vec(end+1) = rho*qtot;
    end
    if j0vec(1) > 0
        j0vec = [0, j0vec];
    end
    
    xx = x0 + j0vec/(rho*qtot) * (maxxqpos - x0);
    j = j0vec(1) + zeros(size(xx));
    
    for i = 1:niter
        j(2:end) = j(1) + rho*intq(xx(2:end));
        j(end) = j0vec(end);
        xx = interp1(j, xx, j0vec);
    end
    
    z = [z, h(xx(2:end-1))];
    x = [x, xx(2:end-1)];
    
end



function [x, z, t] = getTrajectory(x0, z0, u, v, h, dddx, dt)
    % Solves the pde given by
    %   x' = u(x,z)
    %   z' = v(x,z)
    % with initial condition
    %   x(0) = x0
    %   z(0) = z0
    % until z > h
    
    x = x0;
    z = z0;
    t = 0;
    
    while h(x(end)) >= z(end)
        
        xt = x(end);
        zt = z(end);
        
        x(end+1) = xt + u(xt,zt)*dt;
        z(end+1) = zt + v(xt,zt)*dt;
        
        t = t + dt;
        
        if (x(end-1) == x(end)) && (z(end-1) == z(end))
            break
        end
        
    end
    
    if length(x) == 1
        return
    end
    
    a = 0;
    b = dt;
    for i = 1:100
        
        dt = (a+b)/2;
        
        xt = x(end-1);
        zt = z(end-1);
        
        x(end) = xt + u(xt,zt)*dt;
        z(end) = zt + v(xt,zt)*dt;
        
        if h(x(end)) >= z(end)
            a = dt;
        else
            b = dt;
        end
    end
    t = t - dt + b;
end