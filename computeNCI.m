function NCI = computeNCI(estStates,Ps,simulStates)

    %rsme is a vector of size n
    n = size(estStates,1);
    ntimesteps = size(estStates,2);
    
    
    resid = estStates - simulStates;
    Sigm = 0;
    for i = 1:ntimesteps
        Sigm = Sigm + (1/ntimesteps)*resid(:,i)*resid(:,i)';
    end
    smallval = 0.001;
    Sigm = (1-smallval)*Sigm+smallval*eye(n);
    
    NCI = 0;
    for i = 1:ntimesteps
        ep = resid(:,i)'*(Ps(:,:,i)\resid(:,i));
        ep_tilde = resid(:,i)'*(Sigm\resid(:,i));
        NCI = NCI + (10/ntimesteps)*(log10(ep)-log10(ep_tilde));
    end

end