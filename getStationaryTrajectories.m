function [trajectories, times] = getStationaryTrajectories(u, v, h, d, intq, x0, m, kappa, rho, dt, dj, J0, maxxqpos)
    % Computes trajectories for a stationary glacier


    % Construct initial values with uniform flux
    niter = 100;
    [z0vec, x0vec] = equiFlux2(h, d, intq, x0, m, kappa, rho, dj, J0, maxxqpos, niter);
    
    % Compute trajectories
    trajectories = cell(1, length(z0vec));
    times = zeros(1, length(z0vec));
    for i = 1:length(z0vec)
        [x, z, t] = getTrajectory(x0vec(i), z0vec(i), u, v, h, dt);
        trajectories{i} = [x; z];
        times(i) = t;
        plot(x, z, 'k');
    end
    
end




function [z, x] = equiFlux2(h, d, intq, x0, m, kappa, rho, dj, J0, maxxqpos, niter) %(z0vec, u, x0, rho, j0vec, niter)
    
    % Estimate the z-values for x = x0 that gives constant flux between
    % the z's by binary search.
    x = [];
    z = [];
    for j = dj:dj:J0
        za = d(x0);
        zb = h(x0);
        for i = 1:niter
            zi = (za+zb)/2;
            jz = rho*kappa/(m+1) * ((h(0) - d(0)).^(m+1) .* (zi-d(0)) + 1/(m+2) * ((h(0) - zi).^(m+2) - (h(0) - d(0)).^(m+2)));
            if jz < j
                za = zi;
            else
                zb = zi;
            end
        end
        x = [x,x0];
        z = [z,za];
    end
    
    % Estimate the x-values for z = h(x) that gives constant flux between
    % the (x's, z's) by binary search.
    while true
        j = j+dj;
        xa = x(end);
        xb = maxxqpos;
        jxa = rho*intq(xa)+J0;
        jxb = rho*intq(xb)+J0;
        if jxb-jxa < dj
            break
        end
        for i = 1:niter
            xi = (xa+xb)/2;
            jx = rho*intq(xi)+J0;
            if jx < j
                xa = xi;
            else
                xb = xi;
            end
        end
        x = [x,xa];
        z = [z,h(xa)];
    end
    
end



function [x, z, t] = getTrajectory(x0, z0, u, v, h, dt)
    % Solves the pde given by
    %   x' = u(x,z)
    %   z' = v(x,z)
    % with initial condition
    %   x(0) = x0
    %   z(0) = z0
    % until z > h
    % using forward Euler.
    
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