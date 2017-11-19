function [xnew, Pnew] = ekf2TimeUpdateSingle(xold,Pold,Q,del,time,probtype)
    
% description - 
  % does a single time update from time k-1 to time k for the MBR problem 
  % given a model of the system (in the form of a vector/matrix xold, Pold, Q) to 
  % obtain a priori estimates of P and x (P_{k}^- and x_{k}^-) before 
  % measurement information is incorporated.

% input
  % @param xold: n x 1 vector: state vector
  % @param Pold: n x n matrix: covariance matrix
  
  % @param Q: n x n matrix: process noise covariance matrix
  % @param R: j x j matrix: measurement noise covariance matrix
  
  % @param del: scalar: delay between measurement readings.
  % @param time: scalar: time in seconds of time k, used in determining the
  % control action through array alphat in loadParamDataMBR
  
  % @param probtype: string: supports MBR, CSTR, and Bioreactor problem

% output
  % @return xnew: n x 1 vector: a priori state estimate
  % @return Pnew: n x n matrix: a priori covariance matrix estimate
  
  x = xold;
  P = Pold;
  
  n = size(x,1);
  
  switch probtype
      case 'MBR'
        loadParamDataMBR
  
        % numerical integration for one sampling period
        solx = ode45(@(t,x) MBRsimulfun(t,x,alphat(time)),[0,del],x);
        x = deval(solx,del);

        P_col = P(:);
        solP = ode45(@(t,P) MBRsimulCovfun(t,P,alphat(time),Q,n_est,solx),[0,del],P_col); %n_estimate=xold(6)
        P = reshape(deval(solP,del),n,n);
      case 'CSTR'
        % numerical integration for one sampling period
        Q_estimate = 10^5;%10^5;
        f3_estimate = 50;
        solx = ode45(@(t,x) CSTRsimulfun(t,x,Q_estimate,f3_estimate),[0,del],x);
        x = deval(solx,del);

        P_col = P(:);
        solP = ode45(@(t,P) CSTRsimulCovfun(t,P,f3_estimate,Q,solx),[0,del],P_col); %n_estimate=xold(6)
        P = reshape(deval(solP,del),n,n);
      case 'Bioreactor'
      
      otherwise
          warning('invalid probtype')
  end
  
  xnew = x;

  Pnew = P;
  
end