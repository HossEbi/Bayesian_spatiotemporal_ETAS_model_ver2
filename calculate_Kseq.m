% Calculate K
% Mi & Ti  : the vectors of events
% T0=0
% Seismicity forecasting based on a Bayesian spatio?temporal ETAS model
% written by: Hossein Ebrahimian 
% Last update: 11/2022

%% Main

function K = calculate_Kseq (Mi, Ti, Iri, tstart, Ml, theta, Nbo)

intLambda = 0;

No = length(Ti);

for j=2:length(Ti)
    intLambda = intLambda+int3LambdaETAS (Mi, Ti, Iri, Ti(j-1), Ti(j), Ml, [theta(1);1;theta(3:end)]);
end
intLambda = intLambda+int3LambdaETAS (Mi, Ti, Iri, Ti(end), tstart, Ml, [theta(1);1;theta(3:end)]);

K = (No-Nbo)/sum(intLambda);

end

%% END










        




    



















        