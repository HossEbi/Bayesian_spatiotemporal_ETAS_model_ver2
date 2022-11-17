% Generate Mi, Ti for SEQgen (for Ti, using THINNING METHOD)
% Seismicity forecasting based on a Bayesian spatio?temporal ETAS model
% written by: Hossein Ebrahimian and Fatemeh Jalayer  
% Last update: 11/2022

function [Mgen,tgen,numRepeat] = generate_m_time (Mi, Ti, Iri, tstart, tend, Ml, theta, Mmax, mub)

Mgen = generate_M (Ml, theta(1), Mmax);

lambda_max = int2lambdaETAS (Mi,Ti,Iri,tstart,Ml,theta,mub);

accept    = 0;
numRepeat = 0;

tgen = tstart;

while (accept == 0 && tgen <= tend)
    
    numRepeat = numRepeat + 1;
   
    %%% generate time
    tgen = generate_time (tgen, lambda_max);
    
    lambda_gen = int2lambdaETAS (Mi,Ti,Iri,tgen,Ml,theta,mub);
    
    %%% Thinning Algorithm
    ratio = lambda_gen/lambda_max;
    u = rand;
    if (u <= ratio && ratio<=1)
        accept = 1;
    else
        lambda_max = lambda_gen;
    end
    
end

end
