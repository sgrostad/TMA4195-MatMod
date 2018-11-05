%% TestScript:
clear all
close all
clc
%% Accumulation Rate
L = 5;
xf = L;
xs = 2;
N = 101;
q0 = 1;
q = getAccumulationRate(L,xs,xf,N,q0);
figure
hold on
plot(linspace(0,L,N),q,'bo')
plot(xs,0,'*r')
plot(xf,0,'*g')
legend('q','xs','xf')
