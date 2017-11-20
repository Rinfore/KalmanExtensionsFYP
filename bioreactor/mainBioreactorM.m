tic
% number of models
m = 4;

% initialization!
P0 = diag([0.001 0.001 0 0 0 0]); 
R_est = diag([1e-4 1e-4]);
Q_est = diag([0.001 0.001 0 0 0 0]); 
H_mult = [eye(2), zeros(2,4)]; %% a different measurement matrix is required for this multiple-model estimation because it ignores n.

x0_est = [ones(1,4) %x1
    0.01*ones(1,4) %x2
    0.3 0.3 0.12 0.3
    0.25 0.03 0.56 0.25
    0.56 0.55 0.02 0.56
    0.02 0.03 0 0.02];

prior = ones(m,1)/m; %complete uncertainty

% returns a vector of states against time (n by ntimesteps) as the
% first argument. Also returns time series of diagonal of
% covariance matrix, posteriors, states estimated by each model, etc.
[estStatesEKF, EKFP, EKFposterior, estStatesEKFmodels, EKFPmodels, EKFposteriorTimeSeries] = ekf1Multi(prior,x0_est,P0,H_mult,Q_est,R_est,simulMeasur,ntimesteps,del,'Bioreactor',1);
toc

rsmeEKF = computeRSME_Single(estStatesEKF,simulStates);
%rsmeEKF = computeRSME_Multi(estStatesEKF,simulStates(1:end-1,:),estStatesEKFmodels);

toc