clear all
close all
clc

% Created by Riccardo Canola (riccardo.canola@gmail.com)

%% Data

% forgetting factor
lambda = 0.9;

% number of parameters
n = 3;

% samples
m = 1000;

% sampling frequency
fs = 50;

% time vector
t = (0:m-1)/fs;
dt = t(end)/m;

% time-varying parameters [nxm matrix]
for c = 1:n
    
    e = randn;
    f = randn;
    d = randn;
    g = randn;
    
    THETA(c,:) = [linspace(randn,e,m/5),linspace(e,f,m/5),linspace(f,d,m/5),linspace(d,g,m/5),linspace(g,randn,m/5)];  
    % since no flight test data is available I created a random sequence of parameters to test the algorithms
end

% regressors
X = randn(m,n);

% outputs
Z = X*THETA+0.1*randn;

%% Exponentially weighted recursive least squares

[THETA_est_ewrls,COV_ewrls] = EWRLS(X,Z,lambda,n,m);

%% Sequential least squares

f_sls = 1;                                  % parameter estimate update rate (tipically 1 or 2 Hz)
                        
[THETA_est_sls,COV_sls] = SLS(X,Z,lambda,f_sls,t,n,m);
        
%% Frequency domain sequential least squares

f_fsls = 2;                                 % parameter estimate update rate (tipically 1 or 2 Hz)
freq = [0.1:0.04:1.5];                      % frequencies of interest in the discrete Fourier transform

[THETA_est_fsls,COV_fsls] = FSLS(X,THETA,freq,f_fsls,t,dt,n,m);

%% Plots

for i = 1:n
    
    subplot(n,1,i)
    hold on
    grid minor
    plot(t,THETA(i,:),'LineWidth',1);
    plot([t t(end)+dt],THETA_est_ewrls(i,:),'LineWidth',1);
    plot([t t(end)+dt],THETA_est_sls(i,:),'LineWidth',1);
    plot([t t(end)+dt],THETA_est_fsls(i,:),'LineWidth',1);
    legend('true','EWRLS','SLS','FSLS')
    xlabel('$t[s]$','Interpreter','latex','FontSize',18);
    ylabel(sprintf('$\\hat{\\theta}_{%d}$',i),'Interpreter','latex','FontSize',18);
    
end
