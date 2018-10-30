function [EST,COV] = SLS(X,Z,lambda,f,t,n,m)
% SLS - Sequential least squares
%
% INPUTS
% X : matrix of regressors
% Z : outputs(measurements)
% lambda : forgetting factor
% f : parameter estimate update rate [Hz]
% t : vector of time
% n : number of parameters
% m : number of samples
% 
% OUTPUTS
% EST : matrix of estimated parameters
% COV : covariance matrices of parameters estimates
%
% Created by Riccardo Canola (riccardo.canola@gmail.com)

EST = zeros(n,m+1);
count = 1;
s2 = randn;
M = X(1,:)'*X(1,:);
S = X(1,:)'*Z(1,1);
COV(:,:,1) = 10^6 * eye(n,n);

dt= 1/f;                          
m_sls = length(0:dt:t(end));           
w = [1:m/m_sls:m];     

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

    M = M.*lambda+X(k,:)'*X(k,:);
    S = S*lambda+X(k,:)'*Z(k,k);
        
        if ismember(k,w)
            count = count+1;
            EST(:,k+1) = inv(M)*S;
            COV(:,:,count) = s2.*inv(M);
        else
            EST(:,k+1) = EST(:,k);
        end
        
end

end
