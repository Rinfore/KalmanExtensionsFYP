function [x, ernewM]=ekf4CompensMBR(x,ernew,alph,del,time,H)
    
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
  
  er = ernew;

  %n = size(x,1);
  %j = size(y,1);
  
  %H = [eye(4), zeros(4,2)];
      
  loadParamDataMBR
  
  % numerical integration for one sampling period
  solx = ode45(@(t,x) MBRsimulfun(t,er,alphat(time)),[0,del],er);
  er = deval(solx,del);
  
  ernewM = er;
  x = x + alph*H*ernewM;
  
end