% Parameters
[x0, xS, q0, a, m, J0, rho, kappa] = getParam('Engabreen annual');
xF = getStationaryToe(x0, xS, q0, a, J0, rho);

x_initvec = 0.01:0.1:1;
tend = 1;
d = @(x) 0*x;
dddx = @(x) 0;

q    = @(x) getAccumulationRate(x, x0, xS, xF, q0, a);
intq = @(x) getCumulativeAccumulationRate(x, x0, xS, xF, q0, a, J0, rho);
eta = @(x) getStationaryHeightProfile(x, intq, d, m, kappa, J0, rho);
detadx = @(x) getStationaryHeightProfileDerivative(x, q, intq, d, dddx, m, kappa, J0, rho)- dddx(x);

%% charTest
for x_init = x_initvec
charac = clacCharacteristic(x_init, [0, tend], eta, m, kappa);

plot(charac(:,2),charac(:,1))
hold on
end

for tstart = 0.9:-0.1:0.1
    charac = clacCharacteristic(0, [tstart, tend], eta, m, kappa);

plot(charac(:,2),charac(:,1))

end

%% zTest
hold off
x_initvec = 0:0.01:0.95;
k0 = @(x) 300*x*(x-0.2)*exp(-500*(x-0.1)^2)*sin(20*(x-0.1));
kb = @(t) 0;

h = figure;
filename = 'GIFtest.gif';
for tend = 0.01:0.01:0.6
z_vec = [];
x_vec = [];
i=1;
for tstart = tend-0.01:-0.01:0.01
    vals = calcFunctionValueAlongChar(0, [tstart,tend], eta, detadx, m, kappa, k0, kb);
    x_vec(i) = vals(end,1);
    z_vec(i) = vals(end,2);
    i = i+1;
end

for x_init = x_initvec
    vals = calcFunctionValueAlongChar(x_init, [0,tend], eta, detadx, m, kappa, k0, kb);
    x_vec(i) = vals(end,1);
    z_vec(i) = vals(end,2);
    i = i+1;
end
plot(x_vec,z_vec)
axis([0,1.2,-3,3])
drawnow

frame = getframe(h); 
im = frame2im(frame); 
[imind,cm] = rgb2ind(im,256); 
% Write to the GIF File 
if tend==0.01 
  imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.1); 
else 
  imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1); 
end 
end
    
%% bTest (perhaps irrelevant)
hold off
k0 = @(x) 0;
kb = @(t) 0.25*sin(100*t);

h = figure;
filename = 'GIFtest2.gif';

for tend = 0.01:0.05:0.7
z_vec = [];
x_vec = [];
i=1;
for tstart = tend-0.01:-0.01:0.01
    vals = calcFunctionValueAlongChar(0, [tstart,tend], eta, detadx, m, kappa, k0, kb);
    x_vec(i) = vals(end,1);
    z_vec(i) = vals(end,2);
    i = i+1;
end

for x_init = x_initvec
    vals = calcFunctionValueAlongChar(x_init, [0,tend], eta, detadx, m, kappa, k0, kb);
    x_vec(i) = vals(end,1);
    z_vec(i) = vals(end,2);
    i = i+1;
end
plot(x_vec,z_vec)
axis([0,7,-1,1])
drawnow

frame = getframe(h); 
im = frame2im(frame); 
[imind,cm] = rgb2ind(im,256); 
% Write to the GIF File 
if tend==0.01 
  imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
else 
  imwrite(imind,cm,filename,'gif','WriteMode','append'); 
end 
end
    
    
    
    
    