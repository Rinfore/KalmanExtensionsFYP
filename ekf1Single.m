function [xs, Ps] = ekf1Single(x0,P0,H,Q,R,ys,ntimesteps,del,probtype,md,alph,ERCfactor,robustflaglmd)

% description - 
  % performs state estimation for entire ntimesteps for the specified problem given a model of the 
  % system (initial state x0,initial covariance P0 matrix, R matrix, ys matrix of measurements, 
  % and ntimesteps - returning the xs, Ps that contain estimated states and covariance matrices
  % for the entire time horizon)

% input
  
  % @param x0: n x 1 vector: initial state estimate
  % @param P0: n x n matrix: covariance matrix
  % @param H: j x n matrix: output matrix
  % @param Q: n x n matrix: process noise covariance matrix
  % @param R: j x j matrix: measurement noise covariance matrix
  
  % @param ys: m x ntimesteps matrix: measurement vectors comprise the columns of this matrix
  % @param ntimesteps: scalar
  % @param del: scalar: delay between measurement readings.
  
  % @param probtype: string: supports MBR, CSTR, and Bioreactor problem
  % @param alph: j x 1 vector: the value of alpha used for ERC correction
  % @param ERCfactor: scalar: 0 if ERC is not to be used, or the ratio of fast
  % to slow readings
  
  % @param robustflag: binary: 1 if robust measurement update is to be
  % used, 0 otherwise.
  
% output
  % @return xs: n x ntimesteps matrix:  matrix containing diagonal of the a
  % posteriori state estimates, with each column i corresponding to the estimate at time i.
  % @return Ps: n x ntimesteps 3-D matrix: matrix containing diagonal of the a posteriori covariance matrix
  % estimates, with each column i corresponding to the estimate at time i.
  
  n = size(x0,1);
  
  xs = zeros(n,ntimesteps);
  %Ps = zeros(n,ntimesteps);
  Ps = zeros(n,n,ntimesteps);
  
  x = x0;
  P = P0;
  
  for i = 1:ntimesteps
    
    time = ceil(i*del);
    [x, P] = ekf2TimeUpdateSingle(x,P,Q,del,time,probtype,md);
    %x = real(x);
    %P = real(P);
    if ERCfactor
        if mod(i-1,ERCfactor) ~= 0 %%so first data point has a reading
            [x, er] = ekf4CompensMBR(x,er,alph,del,time,H);%%F? resid?
        else
            [x, P, resid, K] = ekf3MeasurUpdate(x,P,R*del,ys(:,i),H,robustflaglmd);
            warning('off','MATLAB:singularMatrix')
            Htransform = (H'*H)\H';
            warning('on','MATLAB:singularMatrix')
            Htransform(isnan(Htransform)) = 0;
            er = (Htransform - K)*resid; %%not sure if this works for other cases
        end
    else
        [x, P, ~, ~] = ekf3MeasurUpdate(x,P,R*del,ys(:,i),H,robustflaglmd);
    end
    %x = real(x);
    %P = real(P);
    xs(:,i) = x;
    Ps(:,:,i) = P;
    %Ps(:,i) = diag(P);
    %Ks(:,i) = K;
    
  end
  
end