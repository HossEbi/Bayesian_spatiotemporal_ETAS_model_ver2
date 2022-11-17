% Generate Ti for SEQgen (THINNING METHOD)
% Seismicity forecasting based on a Bayesian spatio?temporal ETAS model
% written by: Hossein Ebrahimian 
% Last update: 11/2022

function tgen = generate_time (tstart, lambdaMax)

u = rand;

iat = -1/lambdaMax*log(1-u);

tgen = tstart + iat;

end

