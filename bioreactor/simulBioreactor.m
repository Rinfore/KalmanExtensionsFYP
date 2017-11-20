function [simulStates, simulMeasur] = simulBioreactor(Q,R,x0,ntimesteps,del,md)

% description
  % simulates a bioreactor for ntimesteps given noise covariance matrices Q and R
  % and an initial state x0. Returns simulated states for every time step
  % and simulated measurements.

% input
  % @param Q: n x n matrix: (real) process noise covariance matrix
  % @param R: j x j matrix: (real) measurement noise covariance matrix
  % @param x0: n x 1 vector: (real) initial state
  
  % @param ntimesteps: scalar: number of time steps to simulate,
  % non-inclusive of initial state.
  % @param del: scalar: delay between measurement readings.
  
% output
  % @return simulStates: n x ntimesteps matrix: simulated matrix of real states, with each
  % column simulStates_i corresponding to state at time i.
  % @return simulMeasur: j x ntimesteps matrix: simulated matrix of measurements, with each
  % column simulMeasur_i corresponding to measurement at time i.

    x = x0;
    n = size(x,1);
    Q_discrL = Q*del;
    R_discrL = R*del;
    
    pts = 1000;
    ddel = del/pts;

    % % Measure function for the bioreactor problem
    H = [eye(2), zeros(2,4)]; %%hardcoded!
    j = size(H,1); 

    %loadParamData

    simulStates = zeros(n,ntimesteps);
    simulMeasur = zeros(j,ntimesteps);
    
    d_vect = 0.2*ones(1,2000);%[1*10^5*ones(1,500), zeros(1,2000)];%1*10^5*ones(1,2000);%[1*10^5*ones(1,500), zeros(1,2000)];
    cf_vect = 35*ones(1,2000);
    for i = 1:ntimesteps

        %time = ceil(i*del);%%used later to index into the correct control action in the alpha vector
        %Q_discrLpts = Q_discrL/pts;
        %for k = 1:pts
        %    w = mvnrnd(zeros(n,1),Q_discrL);
            
            

        %    dx = [432.9 + 9.68614079478950*10^-2604*x(2)+1.37662*10^8*x(2)*exp(-9057.01/x(1))+4.998*(300-x(1))+34.998*(-x(1)+x(3)) %%correct
        %         4.998*(4-x(2))+34.998*(-x(2)+x(4))-600000*x(2)*exp(-9057.01/x(1))-3000000*x(2)*exp(-6013.95/x(1))
        %         1.37662*10^8*x(4)*exp(-9057.01/x(3)) + 6.49351*10^8*x(4)*exp(-6013.95/x(3)) + 0.001443*1*10^5 + 50*(300 - x(3)) + 13.332*(x(1) - x(3)) 
        %         50/3*(2 - x(4)) + 13.332*(x(2) - x(4)) - 600000*x(4)*exp(-9057.01/x(3)) - 3000000*x(4)*exp(-6013.95/x(3))];
            
            %dx = [432.9 + 9.68614079478950*10^-2604*x(2)+1.37662*10^8*x(2)*exp(-9057.01/x(1))+4.998*(300-x(1))+34.998*(-x(1)+x(3)) %%changed!
            %     4.998*(4-x(2))+34.998*(-x(2)+x(4))-600000*x(2)*exp(-9057.01/x(1))-3000000*x(2)*exp(-6013.95/x(1))
            %     1.37662*10^8*x(4)*exp(-9057.01/x(3)) + 6.49351*10^8*x(4)*exp(-6013.95/x(3)) + 0.001443*1*10^5 + 30*(300 - x(3)) + 13.332*(x(1) - x(3)) 
            %     10*(2 - x(4)) + 13.332*(x(2) - x(4)) - 600000*x(4)*exp(-9057.01/x(3)) - 3000000*x(4)*exp(-6013.95/x(3))];

        %    x = x + dx*ddel+ w'*sqrt(ddel); %%real?????
            %x(1:3) = max(x(1:3),0);
        %end
        
        d_estimate = d_vect(i);
        cf_estimate = cf_vect(i);
        sol = ode45(@(t,x) Bioreactorsimulfun(t,x,d_estimate,cf_estimate,md),[0,del],x);
        x = deval(sol,del);
% 
        w = mvnrnd(zeros(n,1),Q_discrL);
        x = x + w';
        v = mvnrnd(zeros(j,1),R_discrL);
        y = H*x + v';
        %y(1:3) = max(y(1:3),0);

        simulStates(:,i) = x;
        simulMeasur(:,i) = y;

    end
    
end