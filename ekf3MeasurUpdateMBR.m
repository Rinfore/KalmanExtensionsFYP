function [xnewM, PnewM]=ekf3MeasurUpdateMBR(xnew,Pnew,R,y)
    
% description - 
  % does a single measurement update from at time k for the MBR problem 
  % given a model of the system (in the form of a matrix/matrices F, H, 
  % Pnew,xnew,R, y at time k) to obtain a posteriori estimates of P and x 
  % (P_{k}^+ and x_{k}^+) after measurement information is incorporated.
  
% input
  
  % @param xnew: n x 1 matrix: state vector
  % @param Pnew: n x n matrix: covariance matrix
  % @param R: j x j matrix: measurement noise covariance matrix
  
  % @param y: n x 1 vector: measurement
  
% output
  % @return xnewM: n x 1 vector: a posteriori state estimate
  % @return PnewM: n x n matrix: a posteriori covariance matrix estimate
  
  x = xnew;
  P = Pnew;
  n = size(x,1);
  
  H = [eye(4), zeros(4,2)];
      
  resid = y - H*x;
  if ~((sum(sum(isnan(P)))+sum(sum(isinf(P)))) > 0) %ensures no nan and inf in the matrix.
      K = (P*H')/(H*P*H'+R);
      x = x + K*resid;
      P = (eye(n)-K*H)*P*(eye(n)-K*H)'+K*R*K'; % Joseph stabilized version (guarantees positive semidefiniteness)
  end
  %regularization!!!
  %P = (1-1e-4)*P + 1e-4*eye(size(P));%!!!!
  
  xnewM = x;
  PnewM = P;
  
end