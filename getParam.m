function [x0, xS, q0, a, m, J0, rho, kappa] = getParam(type)
    switch type
        case 'Engabreen annual'
            m  = 3.0;
            g  = 9.81;
            rho = 900;
            L = 3400;
            H = 10;
            mu = 1e-20;
            alpha = 15/180*pi;
            Theta = rho*g*H*sin(alpha);
            epsilon = H/L;
            x0 = 0;
            xS = 0.0001;
            xF = 1;
            q0 =  0.1;
            q1 = -1.1;
            a  = (q0-q1)/(xF-xS);
            kappa = (m+2)*(a/2-q0);
            J0 = rho*kappa/(m+2);
            
            year = 3600*24*365;
            Q = 2*epsilon*mu*H*Theta^m/kappa;
            T = H/Q/year
            U = L/T
           
        case 'Engabreen winter'
            m  = 3.0;
            g  = 9.81;
            rho = 900;
            
            L = 3400;
            H = 15;
            mu = 1e-20;
            alpha = 15/180*pi;
            Theta = rho*g*H*sin(alpha);
            epsilon = H/L;
            x0 = 0;
            xS = 0.0001;
            xF = 1;
            q0 =  0.1;
            q1 = -1.1;
            a  = (q0-q1)/(xF-xS);
            kappa = (m+2)*(a/2-q0);
            J0 = rho*kappa/(m+2);
            q0 =  0.5;
            q1 = -0.7;
           
        case 'Engabreen summer'
            m  = 3.0;
            g  = 9.81;
            rho = 900;
            
            L = 3400;
            H = 15;
            mu = 2.4e-24;
            alpha = 15/180*pi;
            Theta = rho*g*H*sin(alpha);
            epsilon = H/L;
            x0 = 0;
            xS = 0.0001;
            xF = 1;
            q0 =  0.1;
            q1 = -1.1;
            a  = (q0-q1)/(xF-xS);
            kappa = (m+2)*(a/2-q0);
            J0 = rho*kappa/(m+2);
            q0 = -0.3;
            q1 = -1.5;
            
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

