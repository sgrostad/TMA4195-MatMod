function [x0, xS, q0, a, m, J0, rho, kappa] = getParam(type)
    switch lower(type)
        case 'engabreen annual'
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
            T = H/Q/year;
            U = L/T;
            fprintf('Time scale T = %.3f years\n', T);
            fprintf('Velocity scale U = %.3f m/year\n', U);
           
        case 'engabreen winter'
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
            
            year = 3600*24*365;
            Q = 2*epsilon*mu*H*Theta^m/kappa;
            T = H/Q/year;
            U = L/T;
            fprintf('Time scale T = %.3f years\n', T);
            fprintf('Velocity scale U = %.3f m/year\n', U);
           
        case 'engabreen summer'
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
            
            year = 3600*24*365;
            Q = 2*epsilon*mu*H*Theta^m/kappa;
            T = H/Q/year;
            U = L/T;
            fprintf('Time scale T = %.3f years\n', T);
            fprintf('Velocity scale U = %.3f m/year\n', U);
            
        otherwise
            warning('Using otherwise in parameters')
            [x0, xS, q0, a, m, J0, rho, kappa] = getParam('engabreen annual');
    end
end

