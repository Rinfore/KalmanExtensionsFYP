function plotprobabilityMulti(EKFposteriorTimeSeries)
    n = size(EKFposteriorTimeSeries,1);
    ntimesteps = size(EKFposteriorTimeSeries,2);
    for i = 1:n
        subplot(n,1,i);
        var1 = plot(1:ntimesteps,EKFposteriorTimeSeries(i,:),'k');
        set(var1,'LineWidth',1);
    end
end