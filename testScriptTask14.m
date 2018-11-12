%% TestScript:
clear all
close all
clc


%% Accumulation Rate

% Parameters
[x0, xS, q0, a, m, J0, rho, kappa] = getParam('Engabreen annual');


% Toe of the glacier
xF = getStationaryToe(x0, xS, q0, a, J0, rho);

% Accumulation rate
q    = @(x) getAccumulationRate(x, x0, xS, xF, q0, a);
intq = @(x) getCumulativeAccumulationRate(x, x0, xS, xF, q0, a, J0, rho);


% Bedrock profile
d    = @(x) 0.0*(1+sin(pi*x/xS + 0));
dddx = @(x) 0.0*pi/xS * cos(pi*x/xS + 0);

% Glacier height profile
h    = @(x) getStationaryHeightProfile(x, intq, d, m, kappa, J0, rho);
dhdx = @(x) getStationaryHeightProfileDerivative(x, q, intq, d, dddx, m, kappa, J0, rho);


% Discretization points
Nx = 1001;
x = linspace(x0,xF+0.1*(xF-x0), Nx);

figure
plot(x, h(x), 'k', 'linewidth', 1)
hold on
plot(x, d(x), 'k', 'linewidth', 1)
axis equal



dt = 0.1;
dj = J0/20;
u = @(x,z) getXVelocity(x, z, kappa, m, h, d);
v = @(x,z) getZVelocity(x, z, kappa, m, h, d, dhdx, dddx);



maxxqpos = getMaxXWithPositiveAccumulation(xS, q0, a);
[c, t] = getStationaryTrajectories(u, v, h, d, dddx, intq, x0, m, kappa, rho, dt, dj, J0, maxxqpos);


for i = 1:length(c)
    xz = c{i}(:,1:5:end);
    %plot(xz(1,:), xz(2,:), 'k');
    %quiver(xz(1,:), xz(2,:), u(xz(1,:), xz(2,:)), v(xz(1,:), xz(2,:)), 0.2, 'k', 'markersize', 200);
end


