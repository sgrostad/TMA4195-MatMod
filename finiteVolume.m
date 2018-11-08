clear all
close all
clc
% Parameters
x0 = 0;
xS = 1;
q0 = 1;
a  = 0.8;
m  = 1.8;
J0 = 1;
rho = 1;
kappa = 1;
[etaInit, xF] = getInitialEta(x0, xS, q0, a, J0, rho, m, kappa);
L = floor(xF * 1.5);
% 
dx = 1e-1;
x = (0:dx:L)';
int = (dx/2:dx:L-dx/2)';
timeMax = 1;
dt = 1e-4;

f = @(eta) kappa/(m+2) * eta.^(m+2);

eta = zeros(length(x), timeMax/dt+1);
eta(:,1) = etaInit(x);

q0 = 0;
xS = 0;
q    = @(x,xF) getAccumulationRate(x, x0, xS, xF, q0, a);
for n = 1:timeMax/dt
    qHat = (q([x(2:end);x(end)+dx],xF) + q(x(1:end),xF)) / 2;
    eta(:,n+1) = eta(:,n) - dt/dx * (f(eta(1:end,n)) - [J0; f(eta(1:end-1,n))]) ...
       + dt*qHat;
    eta(eta(:,n+1)<0,n+1) = 0;
    xF = (find(eta(:,n+1)==0,1) - 1)/length(x)*L;
end