function functionValues = calcFunctionValueAlongChar(x_init, tspan, eta0, deta0dx, m, kappa, k0, kb)


charac = clacCharacteristic(x_init, tspan, eta0, m, kappa);
x = @(t) xFromT(t, charac);

z_init = k0(x_init);
if tspan(1) ~= 0
    z_init = kb(tspan(1));
end
    
[t, k] = ode15s(@(t,k) -kappa*k*(m+1)*eta0(x(t))^m*deta0dx(x(t)), tspan, z_init);

xvec = zeros(length(t),1);
for i = 1:length(t)
    xvec(i) = x(t(i));
end
functionValues = [xvec,k];
end


function x = xFromT(t, charac)

tvec = charac(:,1);
xvec = charac(:,2);

i = 1;
while tvec(i) < t
    i = i+1;
end

if t ~= tvec(i)
    weight1 = (t-tvec(max(1,i-1)))/(tvec(i)-tvec(max(1,i-1)));
    weight2 = (tvec(i)-t)/(tvec(i)-tvec(max(1,i-1)));
    x = weight1*xvec(max(1,i-1))+weight2*xvec(i);
else
    x = xvec(i);
end
end