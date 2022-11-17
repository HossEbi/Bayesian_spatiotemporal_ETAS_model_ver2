% Calculation of rate of events with M=m at t=time using ETAS model (mu_ETAS) for Likelihood Function 
% lambda(t,x,y,M=m|theta,seq,Ml)=beta×exp(-beta(m-Ml))×lambda(t,x,y|theta,seq,Ml)
% Mi & Ti  : the vectors of events
% Seismicity forecasting based on a Bayesian spatio?temporal ETAS model
% written by: Hossein Ebrahimian 
% Last update: 11/2022

function lambda = lambdaETAS_likelihood (Mi, Ti, ri, m, time, Ml, theta, mu_xiyi_Ml)

beta  = theta(1);
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

if ~isempty(index)
    Krt = (q-1)/pi*d.^(2*(q-1))*(p-1)*c^(p-1);
    lambda = beta*exp(-beta*(m-Ml))*(mu_xiyi_Ml+sum((K*exp(alpha*(Mi(index)-Ml))./((time-Ti(index)+c).^p)).*(Krt./(ri{index(end),1}.^2+d.^2).^q)));
else
    lambda = beta*exp(-beta*(m-Ml));
end
    
end






        




    



















        