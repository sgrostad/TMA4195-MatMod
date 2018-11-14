function characteristic = clacCharacteristic(x_init, tspan, eta, m, kappa)

[t,x] = ode15s(@(t,x) kappa * eta(x)^(m+1), tspan, x_init);

characteristic = [t,x];