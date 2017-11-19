function dPdt = MBRMsimulCovfun(t,P,alphatval,Q,n_estimate,sol)

    Area = 1;
    c_10 = 10;
    V0 = 0.1;
    n = size(Q,1);
    
    P = reshape(P, n,n);
    
    x = deval(sol,t); % t is [0,del] (in range of 0 to del inclusive both ends)
    
    A = [2*Area*(1-alphatval)*x(1)*x(3)/(c_10*V0) 0 Area*(1-alphatval)*x(1)^2/(c_10*V0) 0 0 
    -Area*(alphatval)*x(2)*x(3)/(c_10*V0) -Area*(alphatval)*x(1)*x(3)/(c_10*V0) -Area*(alphatval)*x(1)*x(2)/(c_10*V0) 0 0 
    0 0 -Area^(2-n_estimate)*x(3)^(2-n_estimate)*x(5)*(3-n_estimate) 0 -Area^(2-n_estimate)*x(3)^(3-n_estimate) 
    0 0 -Area^(2-n_estimate)*x(3)^(1-n_estimate)*x(5)*(2-n_estimate)*(3-n_estimate)*x(4) -Area^(2-n_estimate)*x(3)^(2-n_estimate)*x(5)*(3-n_estimate) -Area^(2-n_estimate)*x(3)^(2-n_estimate)*(3-n_estimate)*x(4)
    zeros(1,5)]; %%notice that it's zeros(1,5), not zeros(2,6)
    
    dPdt = A*P + P*A' + Q;
    dPdt = dPdt(:);
    
end