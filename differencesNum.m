function etaMat = differencesNum(N,eta0,m,dx,x,dt,kappa,x0,xS,q0,a)
M = length(eta0);
etaMat = zeros(M,N+1);
etaMat(:,1) = eta0;

etaLast = eta0;
for n = 1:N
    %Where is the toe?
    whichxF = find(etaLast==0, 1);
    if isempty(whichxF)
        whichxF = length(etaLast); 
    end
    xF = x(whichxF);
    
    q= (getAccumulationRate(x, x0,xS, xF, q0, a))';
    
    %Forcing first elem. to zero because of BCs
    q(1)=0;
    
    %Computation
    etaPow = etaLast.^(m+1);
    A = dt/dx * kappa *((spdiags(-etaPow,1,M,M))' + spdiags(etaPow,0,M,M));
    A(1,1) = 0;
    etaNew = dt*q + (-A+eye(M))*etaLast;
    
    %eta non-negative:
    etaNew(etaNew<=0) = 0;
    etaMat(:,n+1) = etaNew;
    
    etaLast = etaNew;

  
end

end

