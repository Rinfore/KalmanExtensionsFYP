% MBR

K_est = 2;
K_true = 2;
n_est = 1;
n_true = 1;
Area = 1;
c_10 = 10;
c_20 = 100;
J_0 = 0.06;
%dJ_0 = -4.3200e-04;%n=0
%dJ_0 = -0.0072;%n=1
%dJ_0 = -0.0294;%n=1.5
%dJ_0 = -0.1200;%n=2
V0 = 0.1;

%dJdt = -K_true*Area^(2-n_true)*J_0^(3-n_true);
%d2Jdt2 = -K_true*Area^(2-n_true)*(3-n_true)*J_0^(2-n_true)*dJdt;

alphat = [zeros(1,1) ones(1,20) zeros(1,1) ones(1,5) 20*ones(1,2) zeros(1,20) ones(1,10000) 1];%[ones(1,10000) 1]; %[ones(1,10000) 1];