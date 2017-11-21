% description - 
  % performs MBR simulation for entire ntimesteps of the problem. Please
  % refer to the report for the notation used. e.g. n, j, m, P, R, Q, x, y
  % Allows the calling of mainMBR for the single-model EKF estimation and 
  % mainMBRM for multiple-model EKF estimation.
  % These scripts will return a graph of all estimated parameters in red
  % and real states in blue. rsmeEKF holds the RMSE for each state.

addpath('./cstr')

clear
tic

% duration (in seconds)
dur = 1800/3600; %15 minutes

% time delay between measurements (in seconds)
del = 1/3600;

% number of time steps (minimum)
ntimesteps = floor(dur/del);

% initialize R_real, Q_real, x0_real for MBR simulation
% later
R_real = diag([36 36]);%diag([0.36 0.36]); 
Q_real = diag([13 0.0010 13 0.0010]); %diag([0.36 0.0036 0.36 0.0036]); 
x0_real = [304
    2.5
    303
    2.28];
%[304
%    2.5
%    303
%    2.26];

ERCfactor = false;
alph = NaN;
% simulate the CSTR.
[simulStates, simulMeasur] = simulCSTR(Q_real,R_real,x0_real,ntimesteps,del);
toc

% one is usually commented out. The uncommented one is run and evaluated.

mainCSTR %this is for single-model EKF
%mainCSTRM %this is for multiple-model EKF