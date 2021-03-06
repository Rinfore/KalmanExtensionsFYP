tic
% number of models
m = 4;

% initialization!
P0 = diag([0.001 0.01 0.001 0.001 0.01]); 
R_est = diag([1e-4 1e-4 1e-8 1e-4]);
Q_est = diag([0.001 0.001 0.0001 0.000001 0.001]); 
H_mult = [eye(4), zeros(4,1)]; %% a different measurement matrix is required for this multiple-model estimation because it ignores n.
alph = [diag([5e-5 5e-5 5e-5 5e-5])
      zeros(1,4)];


x0_est = [10
    200
    0.06
    0 %%changed below
    2]; %note that n is not a parameter here.
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


prior = ones(m,1)/m; %complete uncertainty

% returns a vector of states against time (n by ntimesteps) as the
% first argument. Also returns time series of diagonal of
% covariance matrix, posteriors, states estimated by each model, etc.
[estStatesEKF, EKFP, EKFposterior, estStatesEKFmodels, EKFPmodels, EKFposteriorTimeSeries] = ekf1Multi(prior,x0_est,P0,H_mult,Q_est,R_est,simulMeasur,ntimesteps,del,'MBR',0.995,alph,ERCfactor,robustflaglmd,ERCflag);
toc

rsmeEKF = computeRSME_MBR(estStatesEKF,simulStates(1:end-1,:)); %% last row corresponds to n and is irrelevant in this case.
%rsmeEKF = computeRSME_MBRM(estStatesEKF,simulStates(1:end-1,:),estStatesEKFmodels);

NCI = computeNCI(estStatesEKF,EKFP,simulStates(1:end-1,:));

toc