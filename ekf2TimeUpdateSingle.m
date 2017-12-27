function [xnew, Pnew] = ekf2TimeUpdateSingle(xold,Pold,Q,del,time,probtype,md)
    
% description - 
  % does a single time update from time k-1 to time k for either the MBR, 
  % CSTR or bioreactor problem given a model of the system (in the form of
  % a vector/matrix xold, Pold, Q) to obtain a priori estimates of P and x
  % (P_{k}^- and x_{k}^-) before measurement information is incorporated.

% input
  % @param xold: n x 1 vector: state vector
  % @param Pold: n x n matrix: covariance matrix
  
  % @param Q: n x n matrix: process noise covariance matrix

  % @param del: scalar: delay between measurement readings.
  % @param time: scalar: time in seconds of time k, used in determining the
  % control action through array alphat in loadParamDataMBR
  
  % @param probtype: string: supports MBR, CSTR, and Bioreactor problem
  % @param md: scalar: mainly for Bioreactor problem - index of model being
  % used

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
        solP = ode45(@(t,P) CSTRsimulCovfun(t,P,Q,f3_estimate,solx),[0,del],P_col); %n_estimate=xold(6)
        P = reshape(deval(solP,del),n,n);
      case 'Bioreactor'
        % numerical integration for one sampling period
        d_est = 0.2;
        cf_est = 35;
        %md = 4; %%manual override
        
        solx = ode45(@(t,x) Bioreactorsimulfun(t,x,d_est,cf_est,md),[0,del],x);
        x = deval(solx,del);

        P_col = P(:);
        solP = ode45(@(t,P) BioreactorsimulCovfun(t,P,Q,d_est,cf_est,solx,md),[0,del],P_col,odeset('RelTol',1e-2,'AbsTol',1e-4)); %n_estimate=xold(6)
        P = reshape(deval(solP,del),n,n);
        
      otherwise
          warning('invalid probtype')
  end
  
  xnew = x;

  Pnew = P;
  
end