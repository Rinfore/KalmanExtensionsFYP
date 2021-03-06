function dPdt = CSTRsimulCovfun(t,P_col,Q,f3_estimate,sol)

    n = size(Q,1);
    
    P = reshape(P_col, n,n);
    
    x = deval(sol,t); % t is [0,del] (in range of 0 to del inclusive both ends)
    
    %A = [-39.996+(1.24681*10^12*x(2)*exp(-9057.01/x(1)))/x(1)^2, 9.68614079478950*10^-2604+1.37662*10^8*exp(-9057.01/x(1)), 34.998, 0
    %    -((5.43421*10^9*x(2)*exp(-9057.01/x(1)))/x(1)^2)-(1.80419*10^10*x(2)*exp(-6013.95/x(1)))/x(1)^2, -39.996-600000*exp(-9057.01/x(1))-3000000*exp(-6013.95/x(1)), 0, 34.998
    %    13.332, 0, -43.332+(1.24681*10^12*x(4)*exp(-9057.01/x(3)))/x(3)^2+(3.90516*10^12*x(4)*exp(-6013.95/x(3)))/x(3)^2, 1.37662*10^8*exp(-9057.01/x(3))+6.49351*10^8*exp(-6013.95/x(3))
    %    0, 13.332, -((5.43421*10^9*x(4)*exp(-9057.01/x(3)))/x(3)^2)-(1.80419*10^10*x(4)*exp(-6013.95/x(3)))/x(3)^2, -23.332-600000*exp(-9057.01/x(3)) - 3000000*exp(-6013.95/x(3))];
    
    A = [-39.996+(1.24681*10^12*x(2)*exp(-9057.01/x(1)))/x(1)^2+(3.90516*10^12*x(2)*exp(-6013.95/x(1)))/(x(1))^2, 1.37662*10^8*exp(-9057.01/x(1))+6.49351*10^8*exp(-6013.95/x(1)), 34.998, 0
        -((5.43421*10^9*x(2)*exp(-9057.01/x(1)))/x(1)^2)-(1.80419*10^10*x(2)*exp(-6013.95/x(1)))/x(1)^2, -39.996-600000*exp(-9057.01/x(1))-3000000*exp(-6013.95/x(1)), 0, 34.998
        13.332, 0, -13.332-f3_estimate+(1.24681*10^12*x(4)*exp(-9057.01/x(3)))/x(3)^2+(3.90516*10^12*x(4)*exp(-6013.95/x(3)))/x(3)^2, 1.37662*10^8*exp(-9057.01/x(3))+6.49351*10^8*exp(-6013.95/x(3))
        0, 13.332, -((5.43421*10^9*x(4)*exp(-9057.01/x(3)))/x(3)^2)-(1.80419*10^10*x(4)*exp(-6013.95/x(3)))/x(3)^2, -13.332-f3_estimate/3-600000*exp(-9057.01/x(3))-3000000*exp(-6013.95/x(3))];
    
    dPdt = A*P + P*A' + Q;
    dPdt = dPdt(:);
    
end