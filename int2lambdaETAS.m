% Calculation of integral over the whole area A (integral over t,x,y) 
% Int(x,y) lambda(t,x,y|theta,seq,Ml)
% Mi & Ti  : the vectors of events
% Seismicity forecasting based on a Bayesian spatio?temporal ETAS model
% written by: Hossein Ebrahimian 
% Last update: 11/2022

function lambda = int2lambdaETAS (Mi, Ti, Iri, time, Ml, theta, mub)

K     = theta(2);
alpha = theta(3);
c     = theta(4);
p     = theta(5);
d     = theta(6);
q     = theta(7);

index = find(Ti < time);

if length(theta)==8
    gamma = theta(8);
    d = d*exp(gamma*Mi(index));
end

if ~isempty(Iri)
    
    Krt = (q-1)/pi*d.^(2*(q-1))*(p-1)*c^(p-1);
    lambda = mub + sum((K*exp(alpha*(Mi(index)-Ml))./((time-Ti(index)+c).^p)).*(Krt.*Iri(index)));
    
else
    
    Kt = (p-1)*c^(p-1);
    lambda = mub + sum(K*Kt*exp(alpha*(Mi(index)-Ml))./((time-Ti(index)+c).^p));

end

end






        




    



















        