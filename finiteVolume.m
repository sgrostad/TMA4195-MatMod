clear all
close all
clc
% Parameters
x0 = 0;
xS = 1;
L = 2;
q0 = 0.2;
a  = 0.8;
m  = 1.8;
J0 = 11;
rho = 1;
kappa = 1;

dx = 1e-4;
x = 0:dx:L;
timeMax = 1;
dt = 1e-4;
etaInit = zeros(length(x),1);
% etaInit = getInitialEta();
f = @(eta) kappa/(m+2) * eta.^(m+2);

eta = zeros(length(etaInit), timeMax/dt+1);
eta(:,1) = etaInit;

q    = @(x) getAccumulationRate(x, x0, xS, xF, q0, a);
for n = 1:timeMax/dt
    qn = q(x);
    qMid = (qn(2:end) - qn(1:end-1)) / 2;
    qBar = [0; qMid];
    eta(:,n+1) = eta(:,n) - dt/dx * ([0; ...
        f(eta(1:end-1,n)) - f(eta(2:end,n))]) + dt*qBar;
end