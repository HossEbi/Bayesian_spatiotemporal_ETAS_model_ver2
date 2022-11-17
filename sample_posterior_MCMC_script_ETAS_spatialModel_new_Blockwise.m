% This scripts uses MCMC simulation with Metropolis-Hastings procedure to estimate the probability P(Data|Model) 
% Seismicity forecasting based on a Bayesian spatio?temporal ETAS model
% written by: Hossein Ebrahimian and Fatemeh Jalayer  
% Last update: 11/2022

function samples = sample_posterior_MCMC_script_ETAS_spatialModel_new_Blockwise(indexSeq,M,time_MS,tstart,Mc,r,rxy,dA,do_calculate_K,consider_approximate,Nbo,mu_xiyi_Ml,priorPDF_par,THETA)

%% Initialize the Metropolis-Hastings sampler

maxIterations_chain1     = 520;   % Maximum number of iterations
maxIterations_chain2_end = 1000;  % Maximum number of iterations
burnin        = 20;               % burnin period
numChain      = 6;                % number of chain
numUP         = length(THETA);    % number of Uncertain Parameters
kernelType    = 'adaptive';       % kernelType = 'adaptive' / 'nonadaptive'

nbins = 30;                  % number of bins for the histogram

%% Define the PDF for Prior and Proposal Distributions 
   
%%% Prior 
priorPDF = {};
priorPDFfunction = {};
for i=1:numUP 
    priorPDF = [priorPDF,'lognormal'];
    priorPDFfunction = [priorPDFfunction,{@(x,par) lognpdf(x,log(par(1)),par(2))}];  
end
              
%%% Proposal Distributions
proposalPDF = {};
proposalPDF_par = {};
for i=1:numUP 
    proposalPDF = [proposalPDF,'lognormal'];
    proposalPDF_par = [proposalPDF_par,0.30];    % COV
end
weights_prior = [];

%% likelihoodFunction definition

likelihoodFunction = @LikelihoodFunction_spatialModel;

DATA = {M(indexSeq),time_MS(indexSeq),tstart,Mc,r,rxy(indexSeq,:),dA,do_calculate_K,consider_approximate,Nbo,mu_xiyi_Ml};

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Start sampling

nchain = 1;    

while nchain <= numChain
    
    display(['           - Chain Number = ',num2str(nchain)]) 
    
    if nchain==1
        maxIterations = maxIterations_chain1;
    else
        maxIterations = maxIterations_chain2_end;
    end
    
    state  = zeros(numUP,maxIterations);  % Storage space for our samples
    accept = zeros(numUP,maxIterations);  % Storage space for accept/reject decisions

    for iter = 1:maxIterations
        display(['              --- Iteration Number = ',num2str(iter),'/',num2str(maxIterations)]) 
        
        if nchain==1
            for n = 1:numUP
                if (~do_calculate_K || (do_calculate_K && n~=2))
                    [THETA,accept(n,iter)] = sample_posterior_MCMC(THETA,n,proposalPDF{1,n},proposalPDF_par{1,n},priorPDF{1,n},priorPDF_par{1,n},likelihoodFunction,DATA);
                end
            end
            
        else
            
            [THETA,accept(:,iter)] = sample_posterior_MCMC_updated(THETA,seeds,weights,priorPDF,priorPDF_par,weights_prior,likelihoodFunction,DATA);
        
        end

        if do_calculate_K
            if consider_approximate
                Ir = [];
            else
                Ir = calculate_Ir (rxy(indexSeq,:), THETA(6:end), dA, M(indexSeq));
            end
            THETA(2) = calculate_Kseq (M(indexSeq),time_MS(indexSeq), Ir, tstart, Mc, THETA, Nbo);
        end
        
        state(:,iter) = THETA;      
        
    end
        
    %% Samples we take for further analysis
    
    if nchain==1
        seeds = state(:,burnin+1:maxIterations);
    else
        [~,iix] = unique(state','rows');
        seeds = state(:,sort(iix));
        %seeds = state;
    end
    
    if do_calculate_K
        weights = calculateWeights(seeds([1,3:end],:), kernelType);
    else
        weights = calculateWeights(seeds, kernelType);
    end    
 
    nchain = nchain + 1;
    
end

samples = seeds;

%% END


