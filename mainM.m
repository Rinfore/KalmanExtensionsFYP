% description - 
  % performs MBR (or other) simulation for entire ntimesteps of the problem. Please
  % refer to the report for the notation used. e.g. n, j, m, P, R, Q, x, y
  % Allows the calling of mainMBR for the single-model EKF estimation and 
  % mainMBRM for multiple-model EKF estimation.
  % These scripts will return a graph of all estimated parameters in red
  % and real states in blue. rsmeEKF holds the RMSE for each state.

addpath('./mbr')

clear
tic

% duration (in seconds)
dur = 900; %15 minutes

% time delay between measurements (in seconds)
del = 1;

% number of time steps (minimum)
ntimesteps = floor(dur/del);

% initialize R_real, Q_real, x0_real for MBR simulation
% later
R_real = diag([1e-4 1e-4 1e-8 1e-4]); 
Q_real = diag([0.001 0.001 0.0001 0.000001 0 0]); 
x0_real = [10
    200
    0.06
    -0.0072%-0.12%-0.0294%-0.0072
    2
    1];

% simulate the MBR.
[simulStates, simulMeasur] = simulMBR(Q_real,R_real,x0_real,ntimesteps,del);
toc

% one is usually commented out. The uncommented one is run and evaluated.

%mainMBR %this is for single-model EKF
mainMBRM %this is for multiple-model EKF