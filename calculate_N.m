% Calculation of integral of rate of events with M=m in the interval [tstart,tend]
% Seismicity forecasting based on a Bayesian spatio?temporal ETAS model
% written by: Hossein Ebrahimian 
% Last update: 11/2022

% Mi & Ti  : the vectors of events

function N = calculate_N (Mi, Ti, ri, m, tstart, tend, Ml, theta, Int_t_lambda,dA,mu_xy_Ml)

beta  = theta(1);
K     = theta(2);
alpha = theta(3);
c     = theta(4);
p     = theta(5);
d     = theta(6);
q     = theta(7);

index = find(Ti < tend);    

if length(theta)==8
    gamma = theta(8);
    d = d*exp(gamma*Mi(index));
end

Krt = (q-1)/pi*d.^(2*(q-1))*(p-1)*c^(p-1);

if ~isempty(Int_t_lambda) 
    
    if p == 1    
        Io = log((tend-Ti(index)+c)./(Ti(end)-Ti(index)+c));     
    else
        Io = ((tend-Ti(index)+c).^(1-p)-(Ti(end)-Ti(index)+c).^(1-p))/(1-p);
    end
    
    Int_t_lambda_end = sum((K*exp(alpha*(Mi(index)-Ml))).*(Io).*(Krt./(ri(index).^2+d.^2).^q));
    
    N = exp(-beta*(m-Ml))*(mu_xy_Ml*(tend-tstart)+sum(Int_t_lambda)+Int_t_lambda_end)*dA;
    
else
    
    if p == 1
        Io = log((tend-Ti(index)+c)./(tstart-Ti(index)+c));     
    else
        Io = ((tend-Ti(index)+c).^(1-p)-(tstart-Ti(index)+c).^(1-p))/(1-p);
    end

    N = exp(-beta*(m-Ml))*(mu_xy_Ml*(tend-tstart)+sum((K*exp(alpha*(Mi(index)-Ml))).*(Io).*(Krt./(ri(index).^2+d.^2).^q)))*dA;
    
end

end












        




    



















        