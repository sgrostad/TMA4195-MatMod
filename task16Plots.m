clear all
close all
clc
%% plot Volume and Differences together
[x0, xS, q0, a, m, J0, rho, kappa] = getParam('Engabreen winter');
[etaInit, xF, d] = getInitialEta(x0, xS, q0, a, J0, rho, m, kappa);
dx = 1e-2;
dt = 1e-4;
timeMax = 2;
N = timeMax/dt+1;
[x0, xS, q0, a, m, J0, rho, kappa] = getParam('Engabreen annual');
L = max(floor(xF * 1.5),1);
xVol = (dx/2:dx:L-dx/2)';
xDiff = (0:dx:L)';

etaVol = finiteVolume(etaInit, xVol, dx, L, dt, timeMax, xF ,x0, xS, q0, a, J0, m, kappa, rho);
etaDiff = differencesNum(N, etaInit(xDiff), m, dx, xDiff', dt, kappa, x0, xS, q0, a);


figure
createGif = false;
xScale = 3;
if createGif
    filename = 'testAnimated.gif';
    hold off
    p1 = plot(xScale*xVol,etaVol(:,n)+d(xVol),'r');
    set(p1,'LineWidth',2)
    hold on
    p2 = plot(xScale*xDiff,etaDiff(:,n)+d(xDiff),'b');
    set(p2,'LineWidth',1)
    plot(xScale*xVol,d(xVol),'g')
    axis equal
    set(gca,'XTick',[], 'YTick', [])
    set(gca,'FontSize',16)
    legend('Finite Volume','Finite Differences', 'Bedrock')
    gif(filename,'DelayTime',0.1,'LoopCount',5,'frame',gcf);
end

for n = 1:400:length(etaVol)
    hold off
    p1 = plot(xScale*xVol,etaVol(:,n)+d(xVol),'r');
    set(p1,'LineWidth',2)
    hold on
    p2 = plot(xScale*xDiff,etaDiff(:,n)+d(xDiff),'b');
    set(p2,'LineWidth',1)
    plot(xScale*xVol,d(xVol),'g')
    axis equal
    set(gca,'XTick',[], 'YTick', [])
    set(gca,'FontSize',16)
    legend('Finite Volume','Finite Differences', 'Bedrock')
    if createGif
        gif
    else
        pause(0.1)
    end
end

%% Plot varying production
[x0, xS, q0, a, m, J0, rho, kappa] = getParam('Engabreen winter');
[etaInit, xF, d] = getInitialEta(x0, xS, q0, a, J0, rho, m, kappa);
dx = 1e-2;
dt = 1e-4;
timeMax = 3;
L = floor(xF * 1.5);
x = (dx/2:dx:L-dx/2)';
seasonNum = 3;
seasons = {'Engabreen summer','Engabreen annual','Engabreen winter'};
eta = zeros(seasonNum,length(x), floor(timeMax/dt+1));
qSeason = zeros(length(x),seasonNum);
for i = 1:seasonNum
    [x0, xS, q0, a, m, J0, rho, kappa] = getParam(seasons{i});
    [qSeason(:,i), ~] = getAccumulationRate(x, x0, xS, L, q0, a);
    eta(i,:,:) = finiteVolume(etaInit, x, dx, L, dt, timeMax, xF ,x0, xS, q0, a, J0, m, kappa, rho);
    xF = (find(eta(i,:,end)<=1e-15,1) - 1)/length(x)*L;
    if size(xF,2) ~= 1 % No zero elements in eta
        warning('Eta array is full')
        xF = L;
    end
    etaInit = eta(i,:,end);
end

createGif = false;
xScale = 3;
h = figure;
if createGif
    filename = 'testAnimated.gif';
    plot(xScale*x,eta(1,:,1)'+d(x),'r');
    hold on
    plot(xScale*x,d(x),'g')
    plot(xScale*x,qSeason(:,3)/2-1);
    plot(xScale*x,-ones(length(x)),'b--');
    
    xlim([0 xScale*x(end)]);
    ylim([-2.5,2.5]);
    legend('Glacier','Bedrock','q_{season}','q = 0')
    seasonTitle = sprintf('q = %s',seasons{3});
    title(seasonTitle)
    set(gca,'XTick',[], 'YTick', [])
    set(gca,'FontSize',16)
    gif(filename,'DelayTime',0.1,'LoopCount',5,'frame',gcf);
end
for i = 1:seasonNum
    for n = 1:400:size(eta,3)
        hold off
        p1 = plot(xScale*x,eta(i,:,n)'+d(x),'r');
        hold on
        plot(xScale*x,d(x),'g')
        plot(xScale*x,qSeason(:,i)/2-1,'b');
        plot(xScale*x,-ones(length(x)),'b--');
        
        xlim([0 xScale*x(end)]);
        ylim([-2.5,2.5]);
        legend('Glacier','Bedrock','q_{season}','q = 0')
        seasonTitle = sprintf('q = %s',seasons{i});
        title(seasonTitle)
        set(gca,'XTick',[], 'YTick', [])
        set(gca,'FontSize',16)
        drawnow
        if createGif
            gif
        else
            pause(0.1)
        end
    end
end

%% plot increasing glacier
[x0, xS, q0, a, m, J0, rho, kappa] = getParam('Engabreen summer');
[etaInit, xF, d] = getInitialEta(x0, xS, q0, a, J0, rho, m, kappa);
dx = 1e-2;
dt = 1e-4;
timeMax = 1;
N = timeMax/dt+1;
[x0, xS, q0, a, m, J0, rho, kappa] = getParam('Engabreen winter');
L = max(floor(xF * 3),1);
xVol = (dx/2:dx:L-dx/2)';

etaVol = finiteVolume(etaInit, xVol, dx, L, dt, timeMax, xF ,x0, xS, q0, a, J0, m, kappa, rho);


figure
xScale = 3;
for n = 1:400:length(etaVol)
    hold off
    plot(xScale*xVol,etaVol(:,n)+d(xVol),'r');
    hold on
    plot(xScale*xVol,d(xVol),'g')
    legend('Glacier', 'Bedrock')
    xlim([0 xScale*xVol(end)]);
    ylim([-0.5,2.5]);
    set(gca,'XTick',[], 'YTick', [])
    set(gca,'FontSize',16)
    pause(0.3)
end