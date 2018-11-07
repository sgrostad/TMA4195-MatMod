function xF = getStationaryToe(x0, xS, q0, a, J0)

    xF = 1/a * (q0 + a*xS + sqrt(...
        q0^2 + 2*a*(J0/rho + q0*(xS-x0) + a*xS^2 - a*xS)));
    
end