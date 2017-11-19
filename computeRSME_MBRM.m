function rsme = computeRSMEMulti(estStates,simulStates,estStatesmods)
%rsme is a vector of size n
n = size(estStates,1);
ntimesteps = size(estStates,2);
rsme = zeros(n,1);
for i = 1:n %assumes same size for estStates and simulStates
    rsme(i,:) = sqrt(mean((simulStates(i,:) - estStates(i,:)).^2));
end
C = {'k','m','g','c'};

for i = 1:n
    subplot(n,1,i)
    plot(1:ntimesteps,simulStates(i,:),'b')
    hold on
    plot(1:ntimesteps,estStates(i,:),'r')
    for j = 1:size(estStatesmods,2)
        estStatesRelevant = estStatesmods(i,j,:);
        estStatesRelevant = estStatesRelevant(:,:)';
        plot(1:ntimesteps,estStatesRelevant,'color',C{j})
    end
end

end