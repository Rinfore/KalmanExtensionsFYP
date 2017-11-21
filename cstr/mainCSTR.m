% initialization
P0 = diag([0.001 0.001 0.001 0.001]);
R_est = diag([36 36]);
Q_est = diag([13 0.0010 13 0.0010]); 
H = [1 0 0 0
        0 0 1 0];

x0_est = [304
    2.5
    303
    2.28];

% returns a vector of states against time (n by ntimesteps) as the
% first argument.
[estStatesEKF, EKFP] = ekf1Single(x0_est,P0,H,Q_est,R_est,simulMeasur,ntimesteps,del,'CSTR',NaN,alph,ERCfactor);
toc

rsmeEKF = computeRSME_CSTR(estStatesEKF,simulStates);
toc