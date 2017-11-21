function [xs, Ps, posterior, xmodslong, pmodslong, posteriorslong] = ekf1Multi(prior,x0,P0,H,Q,R,ys,ntimesteps,del, probtype, convthreshold, alph, ERCfactor,robustflaglmd)

% description - 
  % performs state estimation for entire ntimesteps for the specified problem given a model of the 
  % system (initial state x0,initial covariance P0 matrix, R matrix, ys matrix of measurements, 
  % and ntimesteps - returning the xs, Ps that contain estimated states and covariance matrices
  % for the entire time horizon). However, dynamics are determined by
  % various models and are are determined by the parameter md passed into
  % the time update. Associated prior probabilities of the model being true
  % are also needed.

% input
  % @param prior: m x 1 vector: vector of probabilities of each model
  % @param x0: n x 1 vector: initial state estimate or
  %          : n x m vector: vector of initial states for each model
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
  m = size(prior,1); %!!!
  
  xs = zeros(n,ntimesteps);
  Ps = zeros(n,ntimesteps);
  
  ignoredIndices = zeros(4,1);
  
  if size(x0,2) == m
      xmods = x0;
  else
      xmods = repmat(x0,[1,m]);
  end
  Pmods = repmat(P0,[1,1,m]);
  Ls = zeros(m,1); %%column vector!!
  xmodslong = zeros(n,m,ntimesteps);
  pmodslong = zeros(n,m,ntimesteps);
  posteriorslong = zeros(m,ntimesteps);
  converged = 0;
  prevLs = ones(m,1);
  if ERCfactor
      ermods = zeros(n,m);
  end
  
  for i = 1:ntimesteps
      
      P_est = reshape(reshape(permute(Pmods,[3 1 2]),m,[]).'*prior,n,[]); % just averaging P's for various models by the probabilities contained in prior
      
      for md = 1:m
         if ~ignoredIndices(md)
             time = ceil(i*del);
             x = xmods(:,md);
             P = Pmods(:,:,md);
             try 
                [x, P] = ekf2TimeUpdateMulti(md,x,P,Q,del,time,probtype); %need to tell which model to follow
                %x = real(x);
                %P = real(P);
                
                if ERCfactor
                    if mod(i-1,ERCfactor) ~= 0 %%so first data point has a reading
                        [x, er] = ekf4CompensMBRM(x,er,alph,del,time,H,md);%%F? resid?
                        Ls(md) = 1; %% arbitrary value that is equal for all models
                    else
                        if converged == 0
                            Ls(md) = computeLikelihoodMulti(H, x, P, del*R, ys(:,i)); %Ls(md) = computeLikelihoodMBRM(H, x, P_est, del*R, ys(:,i));
                        end
                        [x, P, resid, K] = ekf3MeasurUpdate(x,P,R*del,ys(:,i),H,robustflaglmd);
                        warning('off','MATLAB:singularMatrix')
                        Htransform = (H'*H)\H';
                        Htransform(isnan(Htransform)) = 0;
                        warning('on','MATLAB:singularMatrix')
                        er = (Htransform - K)*resid; %%not sure if this works for other cases
                    end
                else
                    if converged == 0
                        Ls(md) = computeLikelihoodMulti(H, x, P, del*R, ys(:,i)); %Ls(md) = computeLikelihoodMBRM(H, x, P_est, del*R, ys(:,i));
                    end
                    
                    [x, P, ~, ~] = ekf3MeasurUpdate(x,P,R*del,ys(:,i),H,robustflaglmd);
                end
                
                %[x, P] = ekf3MeasurUpdate(x,P,del*R,ys(:,i),H);
                %x = real(x);
                %P = real(P);
             catch ME
                 disp(ME)
                 if strcmp(ME.identifier, 'MATLAB:deval:SolOutsideInterval')
                    ignoredIndices(md) = true;
                    Ls(md) = 0;
                    x = zeros(n,1);%or NaN?
                    P = zeros(n,n);%or NaN?
                 else
                    rethrow(ME);
                 end
             end
         else
             Ls(md) = 0;
             x = zeros(n,1);%or NaN?
             P = zeros(n,n);%or NaN?
         end
         
         xmods(:,md) = x;
         Pmods(:,:,md) = P;
         if ERCfactor
             ermods(:,md) = er;
         end
         
         xmodslong(:,md,i) = x;
         pmodslong(:,md,i) = diag(P);
       
      end
      
      %%dealing with NA's
      for md = 1:m
           if (sum(isnan(xmods(:,md)))+sum(isinf(xmods(:,md)))) > 0
               Ls(md) = 0;
               ignoredIndices(md) = true;
           end
      end

      if (ERCfactor && (mod(i-1,ERCfactor)~= 0)) %%times where there are no readings for concentration
          xmodsCleaned = xmods;
          xmodsCleaned(isnan(xmods)) = 0;
          xmodsCleaned(isinf(xmods)) = 0;
          posterior = prior;
          x_est = xmodsCleaned*posterior; % probability-weighted average of state
          P_est = reshape(reshape(permute(Pmods,[3 1 2]),m,[]).'*posterior,n,[]); %%analogous to x_est's calculation.

          xs(:,i) = x_est;
          Ps(:,i) = diag(P_est);
          posteriorslong(:,i) = posterior;
       
          continue
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

     if converged == 0
         if sum(posterior == 0) == m-1 %% if only one model has P = 1, then it's converged.
             converged = 1;
             ignoredIndices = (posterior == 0);
         elseif sum(posterior >= convthreshold) == 1
             converged = 1;
             ignoredIndices = ~(posterior >= convthreshold);
             Ls = prevLs;
         end
     end
     prior = posterior;
      
     %delete low-probability models?
  end
  
end