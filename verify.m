%% Verify the Expression for u
%
%



% Stress tensor
txz  = @(z, htilde) htilde-z;

% Velocity in x-direction
u    = @(z, htilde, d) kappa/(m+1) * (abs(htilde-d).^(m+1) - abs(htilde-z).^(m+1));
dudz = @(z, htilde) kappa*abs(htilde-z).^(m-1).*(htilde-z);

% 


kappa = 3.2;
m = 2.8;
dz = 10.^(-2:-1:-6);



%% Case 1
% d < z < htilde

d = 0.2;
htilde = 0.8;
z = 0.5;
loglog(dz, ...
    (u(z+dz, htilde, d) - u(z-dz, htilde, d))./(2*dz) - ...
    dudz(z, htilde))
hold on;


%% Case 2
% d < htilde < z

d = 0.3;
htilde = 0.5;
z = 0.8;
loglog(dz, ...
    abs((u(z+dz, htilde, d) - u(z-dz, htilde, d))./(2*dz) - ...
    dudz(z, htilde)))


%% Case 3
% htilde < d < z

d = 0.3;
htilde = 0.0;
z = 0.8;
loglog(dz, ...
    abs((u(z+dz, htilde, d) - u(z-dz, htilde, d))./(2*dz) - ...
    dudz(z, htilde)))
hold off


