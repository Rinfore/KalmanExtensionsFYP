function [xs, Ps] = ekf1Single(x0,P0,H,Q,R,ys,ntimesteps,del,probtype)

% description - 
  % performs state estimation for entire ntimesteps for the MBR problem given a model of the 
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
  
% output
  % @return xs: n x ntimesteps matrix:  matrix containing diagonal of the a
  % posteriori state estimates, with each column i corresponding to the estimate at time i.
  % @return Ps: n x ntimesteps 3-D matrix: matrix containing diagonal of the a posteriori covariance matrix
  % estimates, with each column i corresponding to the estimate at time i.
  
  n = size(x0,1);
  
  xs = zeros(n,ntimesteps);
  Ps = zeros(n,ntimesteps);
  
  x = x0;
  P = P0;
  
  for i = 1:ntimesteps
    
    time = ceil(i*del);
    [x, P] = ekf2TimeUpdateSingle(x,P,Q,del,time,probtype);
    %x = real(x);
    %P = real(P);
    [x, P] = ekf3MeasurUpdate(x,P,R*del,ys(:,i),H);
    %x = real(x);
    %P = real(P);
    xs(:,i) = x;
    Ps(:,i) = diag(P);
    %Ks(:,i) = K;
    
  end
  
end