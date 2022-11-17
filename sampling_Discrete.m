%%% Script for sampling from Discrete Function
% Seismicity forecasting based on a Bayesian spatio?temporal ETAS model
% written by: Hossein Ebrahimian 
% Last update: 11/2022

function sampleX = sampling_Discrete(x,px,Nsim)


%% Plot the discrete distribbution

Fx = cumsum(px);

% stairs(x,Fx,'-','color','k','Linewidth',3)
% xlabel('\itX','fontsize',14)
% ylabel('CDF','fontsize',14)
% set(gca,'fontsize',12)
% set(gca,'Xtick',0:10)
% set(gca,'Ytick',0:0.1:1)
% xlim([0 10])
% ylim([0 1])
% grid on

%% Direct Sampling from CDF

sampleX = zeros(1,Nsim);
for i=1:Nsim
    sampleX(i) = x(find(Fx>rand,1,'first'));
end

%% Figure

% figure
% counts = hist(sampleX,x);
% PMF = counts/sum(counts);
% bar(x,PMF,1,'FaceColor',[0.7,0.7,0.8],'EdgeColor',[0.4,0.4,0.4])
% hold on
% bar(x,px,0.50,'k')
% xlabel('\itX','fontsize',14)
% ylabel('PMF','fontsize',14)
% set(gca,'fontsize',12)
% legend(['MCS samples, N =',num2str(Nsim)],'Data')
% 
% figure
% stairs(x,cumsum(PMF),'-','color',[0.7,0.7,0.8],'Linewidth',1.5)
% hold on
% stairs(x,Fx,'-','color','k','Linewidth',3)
% xlabel('\itX','fontsize',14)
% ylabel('CDF','fontsize',14)
% set(gca,'fontsize',12)
% set(gca,'Xtick',0:10)
% set(gca,'Ytick',0:0.1:1)
% xlim([0 10])
% ylim([0 1])
% grid on
% legend(['MCS samples, N =',num2str(Nsim)],'Data')

%%
