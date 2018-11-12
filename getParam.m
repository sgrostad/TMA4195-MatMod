function [x0, xS, q0, a, m, J0, rho, kappa] = getParam(type)
    switch type
        case 'Engabreen annual'
            m  = 3.0;
            g  = 9.81;
            rho = 900;
            
            L = 3400;
            H = 7;
            U = 50/(3600*24*365);
            Q = U*H/L;
            T = L/U;
            mu = 2.4e-24;
            alpha = 15/180*pi;
            Theta = rho*g*H*sin(alpha);
            epsilon = H/L;
            kappa = 2*epsilon*mu*H*Theta^m/Q;
            x0 = 0;
            xS = 0.0001;
            xF = 1;
            q0 =   1 / (Q*3600*24*365*rho/1000);
            q1 = -11 / (Q*3600*24*365*rho/1000);
            a  = (q0-q1)/(xF-xS);
            J0 = rho;
            Q*3600*24*365;
            kappa/(m+2)*H;
           
        case 'Engabreen winter'
            m  = 3.0;
            g  = 9.81;
            rho = 900;
            
            L = 3400;
            H = 7;
            U = 50/(3600*24*365);
            Q = U*H/L;
            T = L/U;
            mu = 2.4e-24;
            alpha = 15/180*pi;
            Theta = rho*g*H*sin(alpha);
            epsilon = H/L;
            kappa = 2*epsilon*mu*H*Theta^m/Q;
            x0 = 0;
            xS = 0.0001;
            xF = 1;
            q0 =  5 / (Q*3600*24*365*rho/1000);
            q1 = -7 / (Q*3600*24*365*rho/1000);
            a  = (q0-q1)/(xF-xS);
            J0 = rho;
            Q*3600*24*365;
            kappa/(m+2)*H;
           
        case 'Engabreen summer'
            m  = 3.0;
            g  = 9.81;
            rho = 900;
            
            L = 3400;
            H = 7;
            U = 50/(3600*24*365);
            Q = U*H/L;
            T = L/U;
            mu = 2.4e-24;
            alpha = 15/180*pi;
            Theta = rho*g*H*sin(alpha);
            epsilon = H/L;
            kappa = 2*epsilon*mu*H*Theta^m/Q;
            x0 = 0;
            xS = 0.0001;
            xF = 1;
            q0 = -3 / (Q*3600*24*365*rho/1000);
            q1 = -15 / (Q*3600*24*365*rho/1000);
            a  = (q0-q1)/(xF-xS);
            J0 = rho;
            Q*3600*24*365;
            kappa/(m+2)*H;
            
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

