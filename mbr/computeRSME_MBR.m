function rsme = computeRSME_MBR(estStates,simulStates)

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
    switch i
        case 1
            axis([0 ntimesteps 0 40])
            title('c1 (mol/m^3)')  
        case 2
            axis([0 ntimesteps -0.3 100 ])
            title('c2 (mol/m^3)')  
        case 3
            axis([0 ntimesteps 0 0.1])
            title('J (mol/m^3)')  
        case 4
            axis([0 ntimesteps -0.05 0.05 ])
            title('dJdt (mol/m^3 t)')  
        case 5
            axis([0 ntimesteps 1 3 ])
            title('K (fouling rate)')  
            if n == 5
                xlabel('time (seconds)')
            end
        case 6
            axis([0 ntimesteps 0 3 ])
            title('n')  
            xlabel('time (seconds)')
        otherwise
            disp('unexpected value of RSME!')
    end
    
end

end