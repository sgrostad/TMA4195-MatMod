function xF = getStationaryToe(x0, xS, q0, a, J0, rho)

    xF = 1/a * (q0 + a*xS + sqrt(...
        q0^2 + 2*a*(J0/rho + q0*(xS-x0))));
%     
%     
%     a0 = a/2;
%     b0 = -q0 - a*xS;
%     c0 = a*xS^2/2 + q0*x0 - J0/rho;
%     
%     xF = 1/(2*a0)* (-b0 + sqrt(b0^2 - 4*a0*c0));
%     
end