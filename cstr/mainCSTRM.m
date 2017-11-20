% number of models
m = 2;

% initialization!
P0 = diag([0.001 0.001 0.001 0.001]);
R_est = diag([36 36]);
Q_est = diag([13 0.0010 13 0.0010]); 
H_mult = [1 0 0 0
        0 0 1 0];

x0_est = [304
    2.5
    303
    2.28];

prior = ones(m,1)/m; %complete uncertainty

% returns a vector of states against time (n by ntimesteps) as the
% first argument. Also returns time series of diagonal of
% covariance matrix, posteriors, states estimated by each model, etc.
[estStatesEKF, EKFP, EKFposterior, estStatesEKFmodels, EKFPmodels, EKFposteriorTimeSeries] = ekf1Multi(prior,x0_est,P0,H_mult,Q_est,R_est,simulMeasur,ntimesteps,del,'CSTR',1);
toc

rsmeEKF = computeRSME_CSTR(estStatesEKF,simulStates); %% last row corresponds to n and is irrelevant in this case.
%rsmeEKF = computeRSMEMulti(estStatesEKF,simulStates(1:end-1,:),estStatesEKFmodels);

toc