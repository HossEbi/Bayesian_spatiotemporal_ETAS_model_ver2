% Function for calculating adaptive/non-adaptive weights in Kernel Density Estimator 
% Seismicity forecasting based on a Bayesian spatio?temporal ETAS model
% written by: Hossein Ebrahimian 
% Last update: 11/2022

function weights = calculateWeights(seeds, type)

[n, nSeeds] = size(seeds);

%% Calculate width (fixed-width) and Md (the number of distinct samples)

% temp = seeds;
% nCopies = 0;
% while size(temp,2)>1
%     temp2 = [];
%     for i = 2:size(temp,2)
%         if isequal(temp(:,i),temp(:,1))
%             nCopies = nCopies+1;
%         else
%             temp2 = [temp2 temp(:,i)];
%         end
%     end
%     temp = temp2;
%     clear temp2;
% end

nCopies = size((unique(seeds','rows'))',2);
if nSeeds==nCopies
    nCopies = 0;
end
    
Md = nSeeds-nCopies;
width = (4/((n+2)*Md))^(1/(n+4));

%% Calculate the local bandwidth factors, lambda (OLD)

% if strcmp(type,'adaptive')
%     alpha = 0.50;   % also [1/n]
%     kp = zeros(nSeeds,1);
%     for i = 1:nSeeds
%         kp(i) = calculateKernel(seeds(:,i),seeds,width*ones(nSeeds,1));
%     end
%     prodkp = prod(kp.^(1/nSeeds));
%     lambda = (prodkp./kp).^alpha;
% else
%     lambda = ones(nSeeds,1);
% end

%% Calculate the local bandwidth factors, lambda (NEW)

if strcmp(type,'adaptive')
    alpha = 0.50;   % also [1/n]
    kp = zeros(nSeeds,1);
    S = cov(seeds');
    display('              ---------- Calculate weights for the samples') 

    for i = 1:nSeeds
        %display(['              ---------- Calculate weight for sample ',num2str(i),'/',num2str(nSeeds)]) 
        X = diag((repmat(seeds(:,i),1,nSeeds)-seeds)'/S*(repmat(seeds(:,i),1,nSeeds)-seeds));
        kp(i) = mean(exp(-0.5*X/width^2)/width^n*1/(sqrt(det(S)*(2*pi)^n)));
    end
    prodkp = prod(kp.^(1/nSeeds));
    lambda = (prodkp./kp).^alpha;
else
    lambda = ones(nSeeds,1);
end

%% Adaptive weights

weights = width*lambda;

return