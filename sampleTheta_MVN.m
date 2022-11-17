%% Function for sampling from a MVN
% Seismicity forecasting based on a Bayesian spatio?temporal ETAS model
% written by: Hossein Ebrahimian 
% Last update: 11/2022

% Mu: should be in a column

function sample_theta = sampleTheta_MVN(Mu,sigma,rho,distribution_type)

n = length(Mu);
beta = diag(sigma);
Cov = beta'*rho*beta;
L = chol(Cov)'; % L = lower triangular matrix L*L'=Cov. Matlab chol function gives the upper triangular matrix L'

sample_theta = Mu + L*randn(n,1);
if strcmp(distribution_type,'lognormal')
      sample_theta = exp(sample_theta);
end

return