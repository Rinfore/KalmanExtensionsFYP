function rsme = computeRSME_Bioreactor(estStates,simulStates)

%rsme is a vector of size n
n = size(estStates,1);
ntimesteps = size(estStates,2);
rsme = zeros(n,1);
for i = 1:n %assumes same size for estStates and simulStates
    rsme(i,:) = sqrt(mean((simulStates(i,:) - estStates(i,:)).^2));
end

f = 2;
for i = 1:f
    subplot(f,1,i)
    var1 = plot(1:ntimesteps,simulStates(i,:),'b');
    set(var1,'LineWidth',1);
    hold on
    var2 = plot(1:ntimesteps,estStates(i,:),'r');
    set(var2,'LineWidth',1);
    switch i
        case 1
            %axis([0 ntimesteps 0 40])
            title('x1 (g/L)')  
        case 2
            %axis([0 ntimesteps -0.3 100 ])
            title('x2 (g/L)')
            xlabel('time (hours)')
        otherwise
            disp('unexpected value of RSME!')
    end
    
end

end