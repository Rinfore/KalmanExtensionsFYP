function [x, ernewM]=ekf4CompensMBR(x,ernew,alph,del,time,H)
    
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
  
% output
  % @return xnewM: n x 1 vector: a posteriori state estimate
  % @return ernewM: n x n matrix: a posteriori residual
  
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