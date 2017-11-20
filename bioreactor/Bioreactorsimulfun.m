function dxdt = Bioreactorsimulfun(t,x,d_est,cf_est,md)

    switch md
        case 1
            dxdt = [(x(3)*x(2)/(x(4)+x(2)) - d_est)*x(1) - x(6)*x(1)
                -x(3)*x(2)/((x(4)+x(2))*x(5))*x(1) + d_est*(cf_est-x(2))
                0
                0
                0
                0];
        case 2
            dxdt = [(x(3)*x(2)/(x(4)*x(1)+x(2)) - d_est)*x(1) - x(6)*x(1)
                -x(3)*x(2)/((x(4)+x(2))*x(5))*x(1) + d_est*(cf_est-x(2))
                0
                0
                0
                0];
        case 3
            dxdt = [(x(3)*x(2) - d_est)*x(1) - x(5)*x(1)
                -x(3)*x(2)/x(2)*x(1) + d_est*(cf_est-x(2))
                0
                0
                0
                0];
        case 4
            dxdt = [(x(3)*x(2)/(x(4)+x(2)) - d_est)*x(1)
                -(x(3)*x(2)/((x(4)+x(2))*x(5))+x(6))*x(1) + d_est*(cf_est-x(2))
                0
                0
                0
                0];
    end

end