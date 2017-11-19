function modelLikelihood = computeLikelihoodMBRM(H, x, P, R, y)

% description -
%   returns the likelihood of the model given the measurement received,
%   covariance matrix, etc. Residuals that are "too far away" are
%   considered low likelihood.

resid = y - H*x;
S = H*P*H' + R;
% making real????
S = real(S);

%regularization!!!
%S = (1-1e-5)*S + 1e-5*eye(size(S));%!!!!
%S = 0.8*S;
%S = 10000*S;
if ((sum(sum(isnan(S))) > 0) || (cond(S) > 10^15))
    modelLikelihood = 0;
else
    modelLikelihood = exp(-resid'*(S\resid/2))/(det(2*pi*S)^(1/2));
    modelLikelihood = real(modelLikelihood);
    if modelLikelihood < 1e-4
        modelLikelihood = 0;
    end
end

end