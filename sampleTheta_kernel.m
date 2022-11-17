%%% Function for sampling from a Kernel
% Seismicity forecasting based on a Bayesian spatio?temporal ETAS model
% written by: Hossein Ebrahimian 
% Last update: 11/2022

%% Main

function sample_theta = sampleTheta_kernel(seeds,weights)

F = cumsum(weights./sum(weights));
%nSeeds = size(seeds,2);
%F = (1:nSeeds)/nSeeds;

%%% Plot
% stairs(1:size(seeds,2),F,'-','color','k','Linewidth',3)
% xlabel('sample','fontsize',14)
% ylabel('CDF','fontsize',14)
% set(gca,'fontsize',12)
% set(gca,'Ytick',0:0.1:1)
% xlim([0 size(seeds,2)])
% ylim([0 1])
% grid on

iKernel = find(F>=rand,1,'first');

Mu = seeds(:,iKernel);
S = (weights(iKernel))^2*cov(seeds');
L = chol(S)';             % L = lower triangular matrix L*L'=S. Matlab chol function gives the upper triangular matrix L'
sample_theta = Mu+L*randn(size(Mu));
 
return

%% END