%% Function for Calculating the Kernel
% Seismicity forecasting based on a Bayesian spatio?temporal ETAS model
% written by: Hossein Ebrahimian 
% Last update: 11/2022

function k = calculateKernel(theta,seeds,weights)

S = cov(seeds');
n = numel(theta);
nSeeds = size(seeds,2);

%% Old Calculation

% k = 0;
% for i = 1:nSeeds
%     k = k + exp(-0.5*(theta-seeds(:,i))'*inv(S)*(theta-seeds(:,i))/weights(i)^2)/...
%         (weights(i)^n*sqrt(det(S)*(2*pi)^n));
% end
% k=k/nSeeds;

%% New Calculation

X = diag((repmat(theta,1,nSeeds)-seeds)'/S*(repmat(theta,1,nSeeds)-seeds));
k = mean(exp(-0.5*X./weights.^2)./(weights.^n)*1/(sqrt(det(S)*(2*pi)^n)));

return