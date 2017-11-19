function dPdt = MBRsimulCovfun(t,P_col,alphatval,Q,n_estimate,sol)

    Area = 1;
    c_10 = 10;
    V0 = 0.1;
    n = size(Q,1);
    
    P = reshape(P_col, n,n);
    
    x = deval(sol,t); % t is [0,del] (in range of 0 to del inclusive both ends)
    
    A = [2*Area*(1-alphatval)*x(1)*x(3)/(c_10*V0) 0 Area*(1-alphatval)*x(1)^2/(c_10*V0) 0 0 0 
    -Area*(alphatval)*x(2)*x(3)/(c_10*V0) -Area*(alphatval)*x(1)*x(3)/(c_10*V0) -Area*(alphatval)*x(1)*x(2)/(c_10*V0) 0 0 0
    0 0 -Area^(2-n_estimate)*x(3)^(2-n_estimate)*x(5)*(3-n_estimate) 0 -Area^(2-n_estimate)*x(3)^(3-n_estimate) -Area^(2-n_estimate)*x(3)^(3-n_estimate)*x(5)*(log(Area) + log(x(3)))
    0 0 -Area^(2-n_estimate)*x(3)^(1-n_estimate)*x(5)*(2-n_estimate)*(3-n_estimate)*x(4) -Area^(2-n_estimate)*x(3)^(2-n_estimate)*x(5)*(3-n_estimate) -Area^(2-n_estimate)*x(3)^(2-n_estimate)*(3-n_estimate)*x(4) Area^(2-n_estimate)*x(3)^(2-n_estimate)*x(5)*x(4)*(1+ (3-n_estimate)*(log(Area) + log(x(3)))) 
    zeros(2,6)];
    
    dPdt = A*P + P*A' + Q;
    dPdt = dPdt(:);
    
end