% description - 
  % performs bioreactor (or other) simulation for entire ntimesteps of the problem. Please
  % refer to the report for the notation used. e.g. n, j, m, P, R, Q, x, y
  % Allows the calling of mainBioreactor for the single-model EKF estimation and 
  % mainBioreactorM for multiple-model EKF estimation.
  % These scripts will return a graph of all estimated parameters in red
  % and real states in blue. rsmeEKF holds the RMSE for each state.

addpath('./bioreactor')

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

% this is the correct R and Q
%R_real = diag([1e-1 1e-1]);
%Q_real = diag([0.01 0.01 0 0 0 0]);

% this is robustness check
R_real = diag([1 1]);
Q_real = diag([1e-2 1e-2 0 0 0 0]);

x0_real = [1 %x1
    0.01 %x2
    0
    0
    0
    0];

md = 3;%1, 2 and 4 work?
switch md
    case 1
        x0_real(3:6,:) = [0.3
            0.25
            0.56
            0.02];
    case 2
        x0_real(3:6,:) = [0.3
            0.03
            0.55
            0.03];
    case 3
        x0_real(3:6,:) = [0.12
            0.56
            0.02
            0];
    case 4    
        x0_real(3:6,:) = [0.3
            0.25
            0.56
            0.02];
    otherwise
        warning('unexpected index of model')
end

ERCflag = false;
ERCfactor = false;
alph = NaN;
robustflaglmd = false;
%robustflaglmd = 2;
% simulate the Bioreactor.
[simulStates, simulMeasur] = simulBioreactor(Q_real,R_real,x0_real,ntimesteps,del,md);
toc

% one is usually commented out. The uncommented one is run and evaluated.

%mainBioreactor %this is for single-model EKF
mainBioreactorM %this is for multiple-model EKF