% New Likelihood Function
% Seismicity forecasting based on a Bayesian spatio?temporal ETAS model
% written by: Hossein Ebrahimian 
% Last update: 11/2022

function [p_teta,lp_teta] = LikelihoodFunction_spatialModel (data, theta)

%% Read parameters

M      = data{1,1};
T      = data{1,2};
tstart = data{1,3};
Ml     = data{1,4};
r      = data{1,5};
rxy    = data{1,6};
dA     = data{1,7};

do_calculate_K        = data{1,8};
consider_approximate  = data{1,9};
Nbo                   = data{1,10};
mu_xiyi_Ml            = data{1,11};

%% Estimate the integral for distance

if consider_approximate
    Ir = [];
else
    Ir =  calculate_Ir (rxy, theta(6:end), dA, M);
end

%% Calculate K 

if do_calculate_K
    theta(2) = calculate_Kseq (M, T, Ir, tstart, Ml, theta, Nbo);
end

%% Main Loop

lprod_lambda = 0;
intlambda = 0;

for i=1:length(T)
    lprod_lambda = lprod_lambda+log(lambdaETAS_likelihood (M, T, r, M(i), T(i), Ml, theta, mu_xiyi_Ml(i)));
    if (i>=2 && ~do_calculate_K)
        intlambda = intlambda+int3LambdaETAS (M, T, Ir, T(i-1), T(i), Ml, theta);
    end
end

if do_calculate_K
    p_teta = exp(lprod_lambda);
    lp_teta = lprod_lambda;
else
    intlambda = intlambda+int3LambdaETAS (M, T, Ir, T(end), tstart, Ml, theta);
    p_teta = exp(lprod_lambda)*exp(-intlambda);
    lp_teta = lprod_lambda-intlambda;
end

end


        




    



















        