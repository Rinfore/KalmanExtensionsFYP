function [xs, Ps, posterior, xmodslong, pmodslong, posteriorslong] = ekf1MBRM(prior,x0,P0,H,Q,R,ys,ntimesteps,del)

% description - 
  % performs state estimation for entire ntimesteps for the MBR problem given a model of the 
  % system (initial state x0,initial covariance P0 matrix, R matrix, ys matrix of measurements, 
  % and ntimesteps - returning the xs, Ps that contain estimated states and covariance matrices
  % for the entire time horizon). However, dynamics are determined by
  % various models and are are determined by the parameter md passed into
  % the time update. Associated prior probabilities of the model being true
  % are also needed.

% input
  % @param prior: m x 1 vector: vector of probabilities of each model
  % @param x0: n x 1 vector: initial state estimate
  % @param P0: n x n matrix: covariance matrix
  % @param H: j x n matrix: output matrix
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
  m = size(prior,1); %!!!
  
  xs = zeros(n,ntimesteps);
  Ps = zeros(n,ntimesteps);
  
  xmods = repmat(x0,[1,m]);
  Pmods = repmat(P0,[1,1,m]);
  Ls = zeros(m,1); %%column vector!!
  xmodslong = zeros(n,m,ntimesteps);
  pmodslong = zeros(n,m,ntimesteps);
  posteriorslong = zeros(m,ntimesteps);
  converged = 0;
  prevLs = ones(m,1);
  
  for i = 1:ntimesteps
      
      P_est = reshape(reshape(permute(Pmods,[3 1 2]),m,[]).'*prior,n,[]); % just averaging P's for various models by the probabilities contained in prior
      
      for md = 1:m
         time = ceil(i*del);
         x = xmods(:,md);
         P = Pmods(:,:,md);
         [x, P] = ekf2TimeUpdateMBRM(md,x,P,Q,del,time); %need to tell which model to follow
         %x = real(x);
         %P = real(P);
         
         if converged == 0
            Ls(md) = computeLikelihoodMBRM(H, x, P, del*R, ys(:,i)); %Ls(md) = computeLikelihoodMBRM(H, x, P_est, del*R, ys(:,i));
         end

         [x, P, ~] = ekf3MeasurUpdateMBRM(x,P,del*R,ys(:,i));
         %x = real(x);
         %P = real(P);
         
         xmods(:,md) = x;
         Pmods(:,:,md) = P;
         
         xmodslong(:,md,i) = x;
         pmodslong(:,md,i) = diag(P);
       
      end
      
       %%dealing with NA's
       for md = 1:m
           if (sum(isnan(xmods(:,md)))+sum(isinf(xmods(:,md)))) > 0
               Ls(md) = 0;
           end
       end
      
      if sum(Ls) == 0
          Ls = prevLs;
      end
      %get the "mean" of x, and P.
      denom = Ls'*prior;
      posterior = Ls.*prior/denom; %column vector element wise with column vector gives column vector
      
      xmodsCleaned = xmods;
      xmodsCleaned(isnan(xmods)) = 0;
      xmodsCleaned(isinf(xmods)) = 0;
      
      x_est = xmodsCleaned*posterior; % probability-weighted average of state
      P_est = reshape(reshape(permute(Pmods,[3 1 2]),m,[]).'*posterior,n,[]); %%analogous to x_est's calculation.
      
      xs(:,i) = x_est;
      Ps(:,i) = diag(P_est);
      posteriorslong(:,i) = posterior;

     if sum(posterior == 0) == m-1 %% if only one model has P = 1, then it's converged.
         converged = 1;
     elseif sum(posterior >= 0.995) == 1
         converged = 1;
         Ls = prevLs;
     end
      prior = posterior;
      
      %delete low-probability models?
  end
  
end