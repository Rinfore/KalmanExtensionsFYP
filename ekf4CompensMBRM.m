function [x, ernewM]=ekf4CompensMBRM(x,ernew,alph,del,time,H,md)
    
% description - 
   % performs estimated residual compensation at time k for the MBR problem 
  
% input
  
  % @param xnew: n x 1 matrix: state vector
  % @param ernew: j x 1 vector: residual
  % @param alph: n x j matrix: alpha matrix (tuning parameter for ERC)
  
  % @param del: scalar: delay between measurement readings.
  % @param time: scalar: time in seconds of time k, used in determining the
  % control action through array alphat in loadParamDataMBR for MBR problem
  
  % @param H: j x n matrix: output matrix
  
  % @param md: scalar: index of the model being used, 1, 2, 3, and 4 
  % corresponding to n = 0, 1, 1.5 and 2 respectively for MBR problem, 
  % n = 1 or 2 for CSTR problem, n = 1, 2, 3, 4 for bioreactor problem.
  
% output
  % @return xnewM: n x 1 vector: a posteriori state estimate
  % @return ernewM: n x n matrix: a posteriori residual
  
  er = ernew;

  %n = size(x,1);
  %j = size(y,1);
  
  %H = [eye(4), zeros(4,2)];
      
  loadParamDataMBR
  switch md
      case 1
        n_estimate = 0;
      case 2
        n_estimate = 1;
      case 3
        n_estimate = 1.5;
      case 4
        n_estimate = 2;
      otherwise
        warning('Unexpected model index.')
  end
  
  % numerical integration for one sampling period
  solx = ode45(@(t,x) MBRMsimulfun(t,er,alphat(time),n_estimate),[0,del],er);
  er = deval(solx,del);
  
  ernewM = er;
  x = x + alph*H*ernewM;
  
end