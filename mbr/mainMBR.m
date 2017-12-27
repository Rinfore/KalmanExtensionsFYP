tic
% initialization
P0 = diag([0.001 0.001 0.001 0.001 0.01 0.01]);
R_est = diag([1e-4 1e-4 1e-8 1e-4]);
Q_est = diag([0.001 0.001 0.0001 0.000001 0.001 0.001]);
H = [eye(4), zeros(4,2)];
alph = [diag([5e-5 5e-5 5e-5 5e-5])
      zeros(2,4)];

x0_est = [10
    200
    0.06
    0 %to be modified later
    2
    1];
x0_est(6) = x0_real(6);
%x0_est(6) = 1;

switch x0_real(6)
    case 0
        x0_est(4) = -4.3200e-04;%n=0
    case 1
        x0_est(4) = -0.0072;%n=0
    case 1.5
        x0_est(4) = -0.0294;%n=0
    case 2.0
        x0_est(4) = -0.1200;%n=0
    otherwise
        warning('unidentified MBR model')
end

% returns a vector of states against time (n by ntimesteps) as the
% first argument.
[estStatesEKF, EKFP] = ekf1Single(x0_est,P0,H,Q_est,R_est,simulMeasur,ntimesteps,del,'MBR',NaN,alph,ERCfactor,robustflaglmd);
toc

rsmeEKF = computeRSME_MBR(estStatesEKF,simulStates);

NCI = computeNCI(estStatesEKF,EKFP,simulStates);

toc