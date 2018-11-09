function [x0, xS, q0, a, m, J0, rho, kappa] = getParam(type)
    switch type
        case 'melting'
            x0 = 0;
            xS = 1;
            q0 = 0;
            a  = 0.8;
            m  = 1.8;
            J0 = 1;
            rho = 1;
            kappa = 1;
        case 'snowing'
            % Change....
            [x0, xS, q0, a, m, J0, rho, kappa] = getParam('melting');
        case 'Sindre history'
            x0 = 0;
            xS = 1;
            q0 = 1;
            a  = 0.8;
            m  = 1.8;
            J0 = 1;
            rho = 1;
            kappa = 1;
        case 'Karine history'
            x0 = 0;
            xS = 0.4;
            q0 = 0.1;
            a  = 2;
            m  = 1.8;
            J0 = 0.26;
            rho = 1;
            kappa = 1;
        otherwise
            [x0, xS, q0, a, m, J0, rho, kappa] = getParam('melting');
    end
end

