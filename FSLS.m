function [EST,COV] = FSLS(X,THETA,freq,f,t,dt,n,m)
% FSLS - Sequential least squares in frequency domain
%
% INPUTS
% X : matrix of regressors
% THETA : matrix of true parameters
% freq : frequencies of interest in the discrete Fourier transform
% f : parameter estimate update rate [Hz]
% t : vector of time
% dt : sampling interval
% n : number of parameters
% m : number of samples
% 
% OUTPUTS
% EST : matrix of estimated parameters
% COV : covariance matrices of parameters estimates
%
% Created by Riccardo Canola (riccardo.canola@gmail.com)

dt_fsls= 1/f;
m_fsls = length(0:dt_fsls:t(end));  
delta_fsls = m/m_fsls;
v = [1:delta_fsls:m];
nfreq = length(freq);

fX = zeros(nfreq,n);
EST = zeros(n,m+1);
COV(:,:,1) = 10^6 * eye(n,n);
count = 1;

noise = 0.1*randn(1,m);
fnoise = fft(noise);        % noise in frequency domain

for k = 0:m-1   
    
    updateomega = 0;

        for q = 0:nfreq-1
            
            omega = (0.1+updateomega)*2*pi;
            fX(q+1,:) = fX(q+1,:)+X(k+1,:)*exp(-j*omega*k*dt);          % recursive Fourier Transform
            updateomega= 0.04*q;
            
        end 
    
    fZ = fX*THETA+fnoise(k+1);      % output in frequency domain
    
        if ismember(k+1,v)
            count = count+1;
            EST(:,k+2) = inv(real(fX'*fX))*real(fX'*fZ(:,k+1));
            s2 = (1/(nfreq-n)).*((fZ(:,k+1)-fX*EST(:,k+2))'*(fZ(:,k+1)-fX*EST(:,k+2)));
            COV(:,:,count) = s2.*inv(real(fX'*fX));
        else
            EST(:,k+2) = EST(:,k+1);
        end
      
end
    
end
