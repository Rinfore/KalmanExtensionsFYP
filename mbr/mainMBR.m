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
    -0.0072
    2
    1];

% returns a vector of states against time (n by ntimesteps) as the
% first argument.
[estStatesEKF, EKFP] = ekf1Single(x0_est,P0,H,Q_est,R_est,simulMeasur,ntimesteps,del,'MBR',NaN,alph,ERCfactor);
toc

rsmeEKF = computeRSME_MBR(estStatesEKF,simulStates);
toc