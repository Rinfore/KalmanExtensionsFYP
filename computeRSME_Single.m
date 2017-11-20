function rsme = computeRSME_Single(estStates,simulStates)

%rsme is a vector of size n
n = size(estStates,1);
ntimesteps = size(estStates,2);
rsme = zeros(n,1);
for i = 1:n %assumes same size for estStates and simulStates
    rsme(i,:) = sqrt(mean((simulStates(i,:) - estStates(i,:)).^2));
end

for i = 1:n
    subplot(n,1,i)
    var1 = plot(1:ntimesteps,simulStates(i,:),'b');
    set(var1,'LineWidth',1);
    hold on
    var2 = plot(1:ntimesteps,estStates(i,:),'r');
    set(var2,'LineWidth',1);
    
    
end

end