%% Function for Calculating the MNV PDF
% Seismicity forecasting based on a Bayesian spatio?temporal ETAS model
% written by: Hossein Ebrahimian 
% Last update: 11/2022

%%% MV Normal: Mu and sigma
%%% MV Lognormal: Mu=median and sigma(logarithmic sigma=COV)

function pdfValue = calculateMVN(theta,Mu,sigma,rho,distribution_type)

n = length(Mu);
beta = diag(sigma);
S = beta'*rho*beta;

if strcmp(distribution_type,'normal') 
    pdfValue = exp(-0.5*(theta-Mu)'/S*(theta-Mu))/sqrt(det(S)*(2*pi)^n);
elseif strcmp(distribution_type,'lognormal')
    lntheta = log(theta);
    lnMu = log(Mu);
    pdfValue = exp(-0.5*(lntheta-lnMu)'/S*(lntheta-lnMu))/(sqrt(det(S)*(2*pi)^n)*prod(theta));
end

return