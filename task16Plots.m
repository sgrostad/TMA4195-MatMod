clear all
close all
clc
%% plot Volume and Differences together
[x0, xS, q0, a, m, J0, rho, kappa] = getParam('Sindre history');
[etaInit, xF, d] = getInitialEta(x0, xS, q0, a, J0, rho, m, kappa);
dx = 1e-2;
dt = 1e-4;
timeMax = 2;
N = timeMax/dt+1;
q0 = 0;
xS = 0;
L = floor(xF * 1.5);
xVol = (dx/2:dx:L-dx/2)';
xDiff = (0:dx:L)';

etaVol = finiteVolume(etaInit, xVol, dx, L, dt, timeMax, xF ,x0, xS, q0, a, J0, m, kappa);
etaDiff = differencesNum(N, etaInit(xDiff), m, dx, xDiff', dt, kappa, x0, xS, q0, a);


figure
for n = 1:400:length(etaVol)
    hold off
    p1 = plot(xVol,etaVol(:,n)+d(xVol),'r');
    set(p1,'LineWidth',2)
    hold on
    axis equal
    p2 = plot(xDiff,etaDiff(:,n)+d(xDiff),'b');
    set(p2,'LineWidth',1)
    plot(xVol,d(xVol),'g')
    legend('Finite Volume','Finite Differences', 'Bedrock')
    pause(0.3)
end

%% Plot varying production
[x0, xS, q0, a, m, J0, rho, kappa] = getParam('Sindre history');
[etaInit, xF, d] = getInitialEta(x0, xS, q0, a, J0, rho, m, kappa);
dx = 1e-2;
dt = 1e-4;
timeMax = 2;
L = floor(xF * 1.5);
x = (dx/2:dx:L-dx/2)';
seasonNum = 5;
seasons = {'melting','Sindre history','melting','snowing','melting'};
eta = zeros(seasonNum,length(x), floor(timeMax/dt+1));
for i = 1:seasonNum
    [x0, xS, q0, a, m, J0, rho, kappa] = getParam(seasons{i});
    eta(i,:,:) = finiteVolume(etaInit, x, dx, L, dt, timeMax, xF ,x0, xS, q0, a, J0, m, kappa);
    etaInit = eta(i,:,end);
end

figure

for i = 1:seasonNum
    for n = 1:400:size(eta,3)
        hold off
        p1 = plot(x,eta(i,:,n)'+d(x),'r');
        hold on
        plot(x,d(x),'g')
        axis equal
        pause(0.05)
    end
end