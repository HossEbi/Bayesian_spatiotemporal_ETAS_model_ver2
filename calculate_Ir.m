% Calculate Ir
% Seismicity forecasting based on a Bayesian spatio?temporal ETAS model
% written by: Hossein Ebrahimian 
% Last update: 11/2022

%% Main

function Ir =  calculate_Ir (rxy, theta, dA, M)

d = theta(1);
q = theta(2);

if length(theta) == 2
    Ir = sum(dA./(rxy.^2+d^2).^q,2);    % dA = dx.dy
else
    gamma = theta(3);
    Ir = sum(dA./(rxy.^2+(repmat(d*exp(gamma*M),1,size(rxy,2))).^2).^q,2);    % dA = dx.dy
end

%% END