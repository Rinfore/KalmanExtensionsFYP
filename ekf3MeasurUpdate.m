function [xnewM, PnewM, resid, K]=ekf3MeasurUpdate(xnew,Pnew,R,y,H,robustflaglmd)
    
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
  % @param H: j x n matrix: output matrix
  % @param robustflag: binary: 1 if robust measurement update is to be
  % used, 0 otherwise.
  
% output
  % @return xnewM: n x 1 vector: a posteriori state estimate
  % @return PnewM: n x n matrix: a posteriori covariance matrix estimate
  % @return resid: j x 1 vector: residual for iteration
  % @return K: n x j matrix: Kalman gain for iteration
  
  x = xnew;
  P = Pnew;
  n = size(x,1);
  j = size(y,1);
      
  resid = y - H*x;
  if ~((sum(sum(isnan(P)))+sum(sum(isinf(P)))) > 0) %ensures no nan and inf in the matrix.
      K = (P*H')/(H*P*H'+R);
      x = x + K*resid;
      if robustflaglmd
          lmd = robustflaglmd;
          gamm = lmd*eigs(inv(inv(P)+H'*(R\H)),1); %%eigs(A,1) returns largest eigenvalue of matrix A
          Re = [R+H*P*H' (P*H')'
              P*H' -gamm^2*eye(n)+P];
          P = (eye(n)-P*[H' eye(n)]*(Re\([H' eye(n)]')))*P;
      else
          P = (eye(n)-K*H)*P*(eye(n)-K*H)'+K*R*K'; % Joseph stabilized version (guarantees positive semidefiniteness)
      end
  else
      disp(P);
  end
  %regularization!!!
  %P = (1-1e-4)*P + 1e-4*eye(size(P));%!!!!
  
  xnewM = x;
  PnewM = P;
  
end