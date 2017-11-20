function dPdt = BioreactorsimulCovfun(t,P_col,Q,d_est,cf_est,sol,md)

    n = size(Q,1);
    
    P = reshape(P_col, n,n);
    
    x = deval(sol,t); % t is [0,del] (in range of 0 to del inclusive both ends)
    
%     switch md
%         case 1
%             A = [-x(6)+(x(3)*x(2))/(-d_est+x(4)+x(2)), -((x(3)*x(1)*x(2))/(-d_est+x(4)+x(2))^2)+(x(3)*x(1))/(-d_est+x(4)+x(2)), (x(1)*x(2))/(-d_est+x(4)+x(2)), -((x(3)*x(1)*x(2))/(-d_est+x(4)+x(2))^2), 0, -x(1)
%                 -((x(3)*x(2))/(x(5)*(x(4)+x(2)))), -d_est+(x(3)*x(1)*x(2))/(x(5)*(x(4)+x(2))^2)-(x(3)*x(1))/(x(5)*(x(4)+x(2))), -((x(1)*x(2))/(x(5)*(x(4)+x(2)))), (x(3)*x(1)*x(2))/(x(5)*(x(4)+x(2))^2), (x(3)*x(1)*x(2))/(x(5)^2*(x(4)+x(2))), 0
%             zeros(4,6)];
%         case 2
%             A = [-x(6)-(x(3)*x(4)*x(1)*x(2))/(-d_est+x(4)*x(1)+x(2))^2+(x(3)*x(2))/(-d_est+x(4)*x(1)+x(2)), -((x(3)*x(1)*x(2))/(-d_est+x(4)*x(1)+x(2))^2)+(x(3)*x(1))/(-d_est+x(4)*x(1)+x(2)), (x(1)*x(2))/(-d_est+x(4)*x(1)+x(2)), -((x(3)*x(1)^2*x(2))/(-d_est+x(4)*x(1)+x(2))^2), 0, -x(1)
%                 -((x(3)*x(2))/(x(5)*(x(4)+x(2)))), -d_est+(x(3)*x(1)*x(2))/(x(5)*(x(4)+x(2))^2)-(x(3)*x(1))/(x(5)*(x(4)+x(2))), -((x(1)*x(2))/(x(5)*(x(4)+x(2)))), (x(3)*x(1)*x(2))/(x(5)*(x(4)+x(2))^2), (x(3)*x(1)*x(2))/(x(5)^2*(x(4)+x(2))), 0
%                 zeros(4,6)];
%         case 3
%             A = [-d_est-x(5)+x(3)*x(2), x(3)*x(1), x(1)*x(2), 0, -x(1), 0
%                 -((x(3)*x(2))/x(4)), -d_est-(x(3)*x(1))/x(4), -((x(1)*x(2))/x(4)), (x(3)*x(1)*x(2))/x(4)^2, 0, 0
%                 zeros(4,6)];
%         case 4
%             A = [(x(3)*x(2))/(-d_est+x(4)+x(2)), -((x(3)*x(1)*x(2))/(-d_est+x(4)+x(2))^2)+(x(3)*x(1))/(-d_est+x(4)+x(2)), (x(1)*x(2))/(-d_est+x(4)+x(2)), -((x(3)*x(1)*x(2))/(-d_est+x(4)+x(2))^2), 0, 0
%                 -x(6)-((x(3)*x(2))/(x(5)*(x(4)+x(2)))), -d_est+x(1)*((x(3)*x(2))/(x(5)*(x(4)+x(2))^2)-x(3)/(x(5)*(x(4)+x(2)))), -((x(1)*x(2))/(x(5)*(x(4)+x(2)))), (x(3)*x(1)*x(2))/(x(5)*(x(4)+x(2))^2), (x(3)*x(1)*x(2))/(x(5)^2*(x(4)+x(2))), -x(1)
%                 zeros(4,6)];
%     end

    switch md
        case 1
            A = [-x(6)+(x(3)*x(2))/(-d_est+x(4)+x(2)), -((x(3)*x(1)*x(2))/(-d_est+x(4)+x(2))^2)+(x(3)*x(1))/(-d_est+x(4)+x(2)), zeros(1,4)
                -((x(3)*x(2))/(x(5)*(x(4)+x(2)))), -d_est+(x(3)*x(1)*x(2))/(x(5)*(x(4)+x(2))^2)-(x(3)*x(1))/(x(5)*(x(4)+x(2))), zeros(1,4)
            zeros(4,6)];
        case 2
            A = [-x(6)-(x(3)*x(4)*x(1)*x(2))/(-d_est+x(4)*x(1)+x(2))^2+(x(3)*x(2))/(-d_est+x(4)*x(1)+x(2)), -((x(3)*x(1)*x(2))/(-d_est+x(4)*x(1)+x(2))^2)+(x(3)*x(1))/(-d_est+x(4)*x(1)+x(2)), zeros(1,4)
                -((x(3)*x(2))/(x(5)*(x(4)+x(2)))), -d_est+(x(3)*x(1)*x(2))/(x(5)*(x(4)+x(2))^2)-(x(3)*x(1))/(x(5)*(x(4)+x(2))), zeros(1,4)
                zeros(4,6)];
        case 3
            A = [-d_est-x(5)+x(3)*x(2), x(3)*x(1), zeros(1,4)
                -((x(3)*x(2))/x(4)), -d_est-(x(3)*x(1))/x(4), zeros(1,4)
                zeros(4,6)];
        case 4
            A = [(x(3)*x(2))/(-d_est+x(4)+x(2)), -((x(3)*x(1)*x(2))/(-d_est+x(4)+x(2))^2)+(x(3)*x(1))/(-d_est+x(4)+x(2)), zeros(1,4)
                -x(6)-((x(3)*x(2))/(x(5)*(x(4)+x(2)))), -d_est+x(1)*((x(3)*x(2))/(x(5)*(x(4)+x(2))^2)-x(3)/(x(5)*(x(4)+x(2)))), zeros(1,4)
                zeros(4,6)];
    end


    dPdt = A*P + P*A' + Q;
    dPdt = dPdt(:);
    
end