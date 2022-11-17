% Calculation of integral in the interval [tstart,tend] and over the whole area A (integral over t,x,y) 
% Int(t) Int(x,y) lambda(t,x,y|theta,seq,Ml)
% Mi & Ti  : the vectors of events
% Seismicity forecasting based on a Bayesian spatio?temporal ETAS model
% written by: Hossein Ebrahimian 
% Last update: 11/2022

%%
function intLambda = int3LambdaETAS (Mi, Ti, Iri, tstart, tend, Ml, theta)

K     = theta(2);
alpha = theta(3);
c     = theta(4);
p     = theta(5);
d     = theta(6);
q     = theta(7);

Kt = (p-1)*c^(p-1);

index = find(Ti < tend); 

if length(theta)==8
    gamma = theta(8);
    d = d*exp(gamma*Mi(index));
end

if p == 1    
    Io = log((tend-Ti(index)+c)./(tstart-Ti(index)+c));     
else
    Io = ((tend-Ti(index)+c).^(1-p)-(tstart-Ti(index)+c).^(1-p))/(1-p);
end

if ~isempty(Iri)
    Kr = (q-1)/pi*d.^(2*(q-1));
    intLambda = sum((K*exp(alpha*(Mi(index)-Ml))).*(Kt*Io).*(Kr.*Iri(index)));
else
    intLambda = sum(K*Kt*exp(alpha*(Mi(index)-Ml)).*Io);
end

end

%% END











        




    



















        