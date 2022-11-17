% Generate Mi for SEQgen
% Seismicity forecasting based on a Bayesian spatio?temporal ETAS model
% written by: Hossein Ebrahimian 
% Last update: 11/2022

function Mgen =  generate_M (Ml, beta, Mmax)

if ~isempty(Mmax)
    FM = 1-exp(-beta*(Mmax-Ml));
else
    FM = 1.0;
end

Mgen = -1/beta*log(1-FM*rand)+Ml;

end
