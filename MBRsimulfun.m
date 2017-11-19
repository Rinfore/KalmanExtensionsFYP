function dxdt = MBRsimulfun(t,x,alphatval)

    Area = 1;
    c_10 = 10;
%    c_20 = 100;
%    J_0 = 0.06;
%    dJ_0 = -0.0072;
    V0 = 0.1;
    
    dxdt = [x(1)^2*Area*x(3)/(c_10*V0)*(1-alphatval) %%changed!
                  -x(1)*x(2)*Area*x(3)/(c_10*V0)*(alphatval)
                  -x(5)*Area^(2-x(6))*x(3)^(3-x(6))
                  -x(5)*Area^(2-x(6))*(3-x(6))*x(3)^(2-x(6))*x(4)
                  0
                  0];
end