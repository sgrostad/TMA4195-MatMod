function etaMat = differencesNum(N,eta0,m,dx,x,dt,kappa,x0,xS0,q0,a)
M = length(eta0);
etaMat = zeros(M,N+1);
etaMat(:,1) = eta0;

%Parameters in q changed permanently:
xS = xS0;
%xS=x0;
%q0=0;


etaLast = eta0;
for n = 1:N
    %Where is the toe?
    whichxF = find(etaLast==0, 1);
    if isempty(whichxF)
        whichxF = length(etaLast); 
    end
    xF = x(whichxF);
    
    %Update q for all positions in time n
    
    %Changing q during the iteration:
    %a_n = min(a+10*dt*n,1);
    %q0_n = max(q0-dt*n, 0);
    %xS_n = max(xS0-dt*n, x0);
    %q= (getAccumulationRate(x, x0,xS_n, xF,q0_n, a))';
    
    %q constant through the iteration (except xF)
    q= (getAccumulationRate(x, x0,xS, xF, q0, a))';
    
    %Forcing first elem. to zero because of BCs
    q(1)=0;
    
    %Computation
    etaPow = etaLast.^(m+1);
    A = dt/dx * kappa *((spdiags(-etaPow,1,M,M))' + spdiags(etaPow,0,M,M));
    A(1,1) = 0;
    etaNew = dt*q + (-A+eye(M))*etaLast;
    
    %eta skal ikke bli negativ:
    etaNew(etaNew<=0) = 0;
    etaMat(:,n+1) = etaNew;
    
    etaLast = etaNew;

  
end


end

