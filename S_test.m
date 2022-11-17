%% This function performs S-Test based on:
% Method 1: Zechar, J. D., Gerstenberger, M. C., & Rhoades, D. A. (2010).
% Likelihood-based tests for evaluating space–rate–magnitude earthquake forecasts.
% Bulletin of the Seismological Society of America, 100(3), 1184-1195.
% Mehod 2: Presented by Ebrahimian and Jalayer (2022)
% Written by: Hossein Ebrahimian
% Last update: 11/2022

%% Main

function [prob1,prob2] = S_test(Nobsxy,Nforexy,Nsimxy)

Nforexy = Nforexy*sum(Nobsxy)/sum(Nforexy);

Likelihood_obs = sum(-Nforexy+Nobsxy.*log(Nforexy)-log(factorial(Nobsxy)));

x = 1:length(Nforexy);

%% Method 1 (Zechar et al. 2010)

simulations = 1000;
Likelihood_sim1 = zeros(simulations,1);
Fx = cumsum(Nforexy/sum(Nforexy));
Nsim = sum(Nobsxy);
for i=1:simulations
    Nxy_sim = zeros(length(Nforexy),1);
    for j=1:Nsim
        sampleX = x(find(Fx>rand,1,'first'));
        Nxy_sim(sampleX) = Nxy_sim(sampleX)+1; 
    end
    Likelihood_sim1(i) = sum(-Nforexy+Nxy_sim.*log(Nforexy)-log(factorial(Nxy_sim)));
end    

%% Method new 

Likelihood_sim2 = zeros(size(Nsimxy,2),1);

for i=1:size(Nsimxy,2)
    Nxy_sim = zeros(length(Nforexy),1);
    Nsimxy(:,i) = Nsimxy(:,i)*sum(Nobsxy)/sum(Nsimxy(:,i));
    Fx = cumsum(Nsimxy(:,i)/sum(Nsimxy(:,i)));
    for j=1:Nsim
        sampleX = x(find(Fx>rand,1,'first'));
        Nxy_sim(sampleX) = Nxy_sim(sampleX)+1; 
    end
    
    Likelihood_sim2(i) = sum(-Nforexy+Nxy_sim.*log(Nforexy)-log(factorial(Nxy_sim)));

end

%% Calculate the probabilities

prob1 = length(find(Likelihood_sim1<=Likelihood_obs))/length(Likelihood_sim1);

prob2 = length(find(Likelihood_sim2<=Likelihood_obs))/length(Likelihood_sim2);

%% END

end

