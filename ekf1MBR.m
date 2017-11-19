function [xs, Ps] = ekf1MBR(x0,P0,Q,R,ys,ntimesteps,del)

% description - 
  % performs state estimation for entire ntimesteps for the MBR problem given a model of the 
  % system (initial state x0,initial covariance P0 matrix, R matrix, ys matrix of measurements, 
  % and ntimesteps - returning the xs, Ps that contain estimated states and covariance matrices
  % for the entire time horizon)

% input
  
  % @param x0: n x 1 vector: initial state estimate
  % @param P0: n x n matrix: covariance matrix
  % @param Q: n x n matrix: process noise covariance matrix
  % @param R: j x j matrix: measurement noise covariance matrix
  
  % @param ys: m x ntimesteps matrix: measurement vectors comprise the columns of this matrix
  % @param ntimesteps: scalar
  % @param del: scalar: delay between measurement readings.
  
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
    [x, P] = ekf2TimeUpdateMBR(x,P,Q,del,time);
    %x = real(x);
    %P = real(P);
    [x, P] = ekf3MeasurUpdateMBR(x,P,R*del,ys(:,i));
    %x = real(x);
    %P = real(P);
    xs(:,i) = x;
    Ps(:,i) = diag(P);
    %Ks(:,i) = K;
    
  end
  
end