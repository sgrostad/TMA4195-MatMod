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
L = max(floor(xF * 1.5),1);
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
qSeason = zeros(length(x),seasonNum);
for i = 1:seasonNum
    [x0, xS, q0, a, m, J0, rho, kappa] = getParam(seasons{i});
    [qSeason(:,i), ~] = getAccumulationRate(x, x0, xS, xF, q0, a);
    eta(i,:,:) = finiteVolume(etaInit, x, dx, L, dt, timeMax, xF ,x0, xS, q0, a, J0, m, kappa);
    xF = (find(eta(i,:,end)<=1e-15,1) - 1)/length(x)*L;
    if size(xF,2) ~= 1 % No zero elements in eta
        warning('Eta array is full')
        xF = L;
    end
    etaInit = eta(i,:,end);
end

createGif = false;
h = figure;
xlim([0 x(end)]);
axis equal
if createGif
    filename = 'testAnimated.gif';
    plot(x,eta(1,:,1)'+d(x),'r');
    hold on
    plot(x,d(x),'g')
    axis equal
    gif(filename,'DelayTime',0.05,'LoopCount',5,'frame',gcf);
end
for i = 1:seasonNum
    for n = 1:400:size(eta,3)
        hold off
        p1 = plot(x,eta(i,:,n)'+d(x),'r');
        hold on
        plot(x,d(x),'g')
        plot(x,qSeason(:,i),'b');
        axis equal
        drawnow
        if createGif
            gif
        end
    end
end