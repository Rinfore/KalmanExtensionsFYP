tic
% initialization
P0 = diag([0.001 0.001 0 0 0 0]);
R_est = diag([1e-4 1e-4]);
Q_est = diag([0.001 0.001 0 0 0 0]);
H = [eye(2), zeros(2,4)];

x0_est = [1 %x1
    0.01 %x2
    0
    0
    0
    0];

%md = 4; %%only for manual override
switch md
    case 1
        x0_est(3:6,:) = [0.3
            0.25
            0.56
            0.02];
    case 2
        x0_est(3:6,:) = [0.3
            0.03
            0.55
            0.03];
    case 3
        x0_est(3:6,:) = [0.12
            0.56
            0.02
            0];
    case 4    
        x0_est(3:6,:) = [0.3
            0.25
            0.56
            0.02];
    otherwise
        warning('unexpected index of model')
end

% returns a vector of states against time (n by ntimesteps) as the
% first argument.
[estStatesEKF, EKFP] = ekf1Single(x0_est,P0,H,Q_est,R_est,simulMeasur,ntimesteps,del,'Bioreactor',md,alph,ERCfactor,robustflaglmd);
toc

rsmeEKF = computeRSME_Single(estStatesEKF,simulStates);
toc