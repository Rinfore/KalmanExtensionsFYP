function [simulStates, simulMeasur] = simulMBR(Q,R,x0,ntimesteps,del,ERCfactor,factorInd,EulerMaru)

% description
  % simulates an MBR for ntimesteps given noise covariance matrices Q and R
  % and an initial state x0. Returns simulated states for every time step
  % and simulated measurements.

% input
  % @param Q: n x n matrix: (real) process noise covariance matrix
  % @param R: j x j matrix: (real) measurement noise covariance matrix
  % @param x0: n x 1 vector: (real) initial state
  
  % @param ntimesteps: scalar: number of time steps to simulate,
  % non-inclusive of initial state.
  % @param del: scalar: delay between measurement readings.
  % @param EulerMaru: binary: 1 if euler-maruyama is to be used
  
% output
  % @return simulStates: n x ntimesteps matrix: simulated matrix of real states, with each
  % column simulStates_i corresponding to state at time i.
  % @return simulMeasur: j x ntimesteps matrix: simulated matrix of measurements, with each
  % column simulMeasur_i corresponding to measurement at time i.

    x = x0;
    n = size(x,1);
    Q_discrL = Q*del;
    R_discrL = R*del;
    
    pts = 500;
    ddel = del/pts;

    % % Measure function for the MBR problem
    H = [eye(4), zeros(4,2)]; %%hardcoded!
    j = size(H,1); 

    loadParamData

    simulStates = zeros(n,ntimesteps);
    simulMeasur = zeros(j,ntimesteps);

    for i = 1:ntimesteps

        time = ceil(i*del);%%used later to index into the correct control action in the alpha vector
        %%%%Q_discrLpts = Q_discrL/pts;
        if EulerMaru
             for k = 1:pts
                 w = mvnrnd(zeros(n,1),Q_discrL);


                 dx = [x(1)^2*Area*x(3)/(c_10*V0)*(1-alphat(time)) %%changed!
                       -x(1)*x(2)*Area*x(3)/(c_10*V0)*(alphat(time))
                       -x(5)*Area^(2-x(6))*x(3)^(3-x(6))
                       -x(5)*Area^(2-x(6))*(3-x(6))*x(3)^(2-x(6))*x(4)
                       0
                       0];

                 x = x + dx*ddel+ w'*sqrt(ddel);

                 x(1:3) = max(x(1:3),0);

             end
        else
             sol = ode45(@(t,x) MBRsimulfun(t,x,alphat(time)),[0,del],x);
             x = deval(sol,del);

             w = mvnrnd(zeros(n,1),Q_discrL);
             x = x + w';
        end
        v = mvnrnd(zeros(j,1),R_discrL);
        y = H*x + v';
        y(1:3) = max(y(1:3),0);
        
        %%begin added%%
        % special code for ERC
        if ERCfactor %%if ERCfactor is 0, runs with regular sampling
            if mod(i-1,ERCfactor) ~= 0 %%so first data point has a reading
                y(factorInd) = NaN;
            end
        end
        %%end added%%

        simulStates(:,i) = x;
        simulMeasur(:,i) = y;

    end
    
end