%% TestScript:
clear all
close all
clc


%% Accumulation Rate

% Parameters
x0 = 0;
xS = 1;
q0 = 0.2;
a  = 0.8;
m  = 1.8;
J0 = 11;
rho = 1;
kappa = 1;


% Toe of the glacier
[~, xF] = getAccumulationRate(0, x0, xS, q0, a, J0, rho);

% Accumulation rate
q    = @(x) getAccumulationRate(x, x0, xS, q0, a, J0, rho);
intq = @(x) getCumulativeAccumulationRate(x, x0, xS, q0, a, J0, rho);

% Bedrock profile
d    = @(x) 0.3*(1-sin(2*pi*x/(xF-x0)));
dddx = @(x) 0.1*2*pi/(xF-x0) * cos(2*pi*x/(xF-x0));

% Glacier height profile
h    = @(x) getStationaryHeightProfile(x, q, intq, d, dddx, m, kappa, J0, rho);
dhdx = @(x) getStationaryHeightProfileDerivative(x, q, intq, d, dddx, m, kappa, J0, rho);


% Discretization points
Nx = 1001;
x = linspace(x0,xF+0.1*(xF-x0), Nx);

figure
plot(x,h(x))
hold on
plot(x,d(x))
axis equal


% 
dt = 0.1;
dj = 0.1;
u = @(x,z) getXVelocity(x, z, kappa, m, h, d);
v = @(x,z) getZVelocity(x, z, kappa, m, h, d, dhdx, dddx);
c = getStationaryTrajectories(u, v, h, d, x0, m, kappa, rho, dt, dj);

for i = 1:size(c,1)
    x = c{1};
    plot(x(1,:), x(2,:));
end
