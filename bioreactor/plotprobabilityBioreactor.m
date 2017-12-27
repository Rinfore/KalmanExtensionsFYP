function plotprobabilityBioreactor(EKFposteriorTimeSeries)
    n = size(EKFposteriorTimeSeries,1);
    ntimesteps = size(EKFposteriorTimeSeries,2);
    for i = 1:n
        subplot(n,1,i);
        var1 = plot(1:ntimesteps,EKFposteriorTimeSeries(i,:),'k');
        set(var1,'LineWidth',1);
        switch i
            case 1
                title('Probability of Model 1')  
            case 2
                title('Probability of Model 2')
            case 3
                title('Probability of Model 3') 
            case 4
                title('Probability of Model 4') 
            otherwise
                disp('unexpected value of i!')
        end
    end
end