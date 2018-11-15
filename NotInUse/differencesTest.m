%Test of finite differences solver.

clear all
%close all
%clc


%% Accumulation Rate

% Parameters
x0 = 0;
xS = 0.4;
q0 = 0.1;
a  = 2;
m  = 1.8;
J0 = 0.26;
rho = 1;
kappa = 1;

%Using steady state as initial state
% Toe of the glacier
xF = getStationaryToe(x0, xS, q0, a, J0, rho);

% Accumulation rate
q    = @(x) getAccumulationRate(x, x0, xS, xF, q0, a);
intq = @(x) getCumulativeAccumulationRate(x, x0, xS, xF, q0, a, J0, rho);

% Bedrock profile
%d    = @(x) 0.1*(1-sin(2*pi*x/(xF-x0)));
%dddx = @(x) 0.1*2*pi/(xF-x0) * cos(2*pi*x/(xF-x0));
d = @(x) 0;
dddx = @(x) 0;

% Glacier height profile
h    = @(x) getStationaryHeightProfile(x, intq, d, m, kappa, J0, rho);
dhdx = @(x) getStationaryHeightProfileDerivative(x, q, intq, d, dddx, m, kappa, J0, rho);

% Discretization points i x-direction
M = 500;
x = linspace(x0,1, M);
dx = x(2)-x(1);

%Discretization in time
N = 2000;
tend = 2;
dt = tend/N;
%initial values: 
eta0 = (h(x)-d(x))';

%make heights over time
etaMat = differencesNum(N,eta0,m,dx,x,dt,kappa,x0,xS,q0,a);
%bottom=d(x);
%Plot
for i=1:100:N+1
    figure(1)
    plot(x, etaMat(:,i))
    hold on;
    %drawnow()
    %pause(0.2)
end
