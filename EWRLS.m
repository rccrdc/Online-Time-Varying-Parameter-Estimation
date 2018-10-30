function [EST,COV] = EWRLS(X,Z,lambda,n,m)
% EWRLS - Exponentially weighted least squares
%
% INPUTS
% X : matrix of regressors
% Z : outputs(measurements)
% lambda : forgetting factor
% n : number of parameters
% m : number of samples
% 
% OUTPUTS
% EST : matrix parameters estimates
% COV : covariance matrices of estimates
%
% Created by Riccardo Canola (riccardo.canola@gmail.com)

COV(:,:,1) = 10^6 * eye(n,n);
s2 = randn;
D(:,:,1) = 10^6 * eye(n,n);
EST = zeros(n,m+1);

for k = 1:m

    v = Z(k,k)-X(k,:)*EST(:,k);

        if k <= 5*n
            s2 = 1/(k+1)*(k*s2+v^2);
        end

        if k == 5*n
            s2 = 1/(k+1-n)*(k*s2+v^2);
        end

        if k > 5*n
            s2 = 1/(k+1-n)*((k-n)*s2+v^2);
        end

    K = D(:,:,k)*X(k,:)'./(lambda+X(k,:)*D(:,:,k)*X(k,:)');
    D(:,:,k+1) = (D(:,:,k)-K*(X(k,:)*D(:,:,k))).*(1/lambda);
    EST(:,k+1) = EST(:,k)+K*(Z(k,k)-X(k,:)*EST(:,k));

    COV(:,:,k+1) = s2.*D(:,:,k+1);

end

end

