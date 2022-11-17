% Seismicity forecasting based on a Bayesian spatio?temporal ETAS model
% written by: Hossein Ebrahimian and Fatemeh Jalayer  
% Last update: 11/2022

% -----------------------------------------------------
% If you use this code, please cite the fillowing two articles:
% (1) Ebrahimian, H., Jalayer, F., Maleki Asayesh, B., Hainzl, S., Zafarani, H. Improvements to seismicity forecasting based on a Bayesian spatio?temporal ETAS model. Scientific Reports, (2022) https://doi.org/10.1038/s41598-022-24080-1.
% (2) Ebrahimian, H. & Jalayer, F. Robust seismicity forecasting based on Bayesian parameter estimation for epidemiological spatio-temporal aftershock clustering models. Scientific Reports 7, 9803 (2017), https://doi.org/10.1038/s41598-017-09962-z.
% -----------------------------------------------------

% Note: This code is derived based on the definitions that are presented in (1) above. Therefore, for some input parameters, it is required to see the paper accordingly.  

clc
clear

%% Important Input Parameters

%%% Forecasting interval (day)
T_start = 11.6250;
T_end   = 12.00;

%%% Method of Analysis
%method_for_analysis = 'Fast';
%method_for_analysis = 'Semi-Fast';
method_for_analysis = 'Slow';

%%% Catalog
filename = 'catalog.txt';  

%%% Magnitude
Mc   = 3.4;     % Lower-bound Magnitude for drawing the number of events within the Catalog
Mmax = 7.5;     % Upper-bound Magnitude (we do not define any upper bound) - you can define any if you want
vecM = [Mc,4,5,6,7];  

%%% Output Dir
output_Dir = 'output';

%%% Control Parameter [0 or 1]
use_background          = 1;    % (=1) consider background seismicity that is in background seismicity.txt file; (=0) consider BGS=0  
use_kernel_distance_m   = 1;    % (=1) magnitude-dependent spatial kernel density; (=0) simple spatial kernel density

%%% Aftershock Zone
deltaGrid = 0.01;
lonMin = 45;
lonMax = 47;
latMin = 32.50;
latMax = 35.50;

%%% ETAS initial parameters [beta,K,alfa,c,p,d,q]
beta_ini  = 1.0*log(10);   
K_ini     = 5.0;          % if do_calculate_K = 1, you can assign an arbitrary value (it is not considered)
alpha_ini = 1.0*log(10);   
c_ini     = 10^(-1.53);
p_ini     = 1.10;
d_ini     = 1.00;
q_ini     = 1.50;
gamma_ini = 0.20;

%%% COV
cov_beta  = 0.5;
cov_K     = 1.0;          % if do_calculate_K = 1, you can assign an arbitrary value (it is not considered)
cov_alpha = 0.5;
cov_c     = 0.5;     
cov_p     = 0.5;
cov_d     = 0.5;
cov_q     = 0.5;
cov_gamma = 0.5;

%%% Vector of ETAS parameters
vec_beta  = 0.01:0.01:4.00;
vec_K     = 0.1:0.1:20;
vec_alpha = 0.01:0.01:4.00;
vec_c     = 0.001:0.001:0.50; 
vec_p     = 0.01:0.01:4.00;
vec_d     = 0.01:0.01:15.00;
vec_q     = 0.01:0.01:5.00;
vec_gamma = 0.01:0.01:1.00;

%% Control Pannel

do_MCMC_updating             = 0;   

do_find_Robust_estimate      = 0;

do_N_test_S_test             = 1;

%% MAIN ******************************************************************

%% Method of Analysis details

switch method_for_analysis
    
    case 'Fast'
        do_calculate_K = 0;        % =1 Calculate K / =0 Leave it to MCMC 
        consider_approximate = 1;  % =1 Approximate / =0 Exact for estimating Int(x,y)   
    case 'Semi-Fast' 
        do_calculate_K = 1;        % =1 Calculate K / =0 Leave it to MCMC 
        consider_approximate = 1;  % =1 Approximate / =0 Exact for estimating Int(x,y)   
    case 'Slow'
        do_calculate_K = 1;        % =1 Calculate K / =0 Leave it to MCMC 
        consider_approximate = 0;  % =1 Approximate / =0 Exact for estimating Int(x,y) 
    otherwise 
        warning('No method for analysis is defined!')
        
end

%% Initial Assignments for MCMC (the order of parameters is important)

priorPDF_par = {[beta_ini,cov_beta],...
                [K_ini,cov_K],...
                [alpha_ini,cov_alpha],...
                [c_ini,cov_c],...
                [p_ini,cov_p],...
                [d_ini,cov_d],...
                [q_ini,cov_q]};
    
THETA = [beta_ini;K_ini;alpha_ini;c_ini;p_ini;d_ini;q_ini];  % Initial Start of THETA      

if use_kernel_distance_m == 1
    priorPDF_par = [priorPDF_par,[gamma_ini,cov_gamma]];
    THETA = [THETA;gamma_ini];
end

%% Read The Catalog

Catalog = load(filename);

latitude  = Catalog(:,1);
longitude = Catalog(:,2);
M         = Catalog(:,3);
time      = Catalog(:,4);
time_MS   = Catalog(:,5);

tstart = T_start-time(1);
tend   = T_end-time(1); 

%% Define the grid and find the distance of the center of each cell unit from events in seq

Xgrid = lonMin:deltaGrid:lonMax;
Ygrid = latMin:deltaGrid:latMax;

Xcgrid = Xgrid(1:end-1)+deltaGrid/2;
Ycgrid = Ygrid(1:end-1)+deltaGrid/2;

num_grid = length(Xcgrid)*length(Ycgrid);

Ggrid = zeros(num_grid,2);
count=1;
for j=1:length(Ycgrid)
    for k=1:length(Xcgrid)
        Ggrid(count,:) = [Ycgrid(j),Xcgrid(k)];
        count=count+1;
    end
end
rxy = calculate_rxy(latitude,longitude,Ggrid);   % rxy = [length(Seq),length(Ggrid)];

dA = prod(topgeo(Ygrid(2)*pi/180,Xgrid(2)*pi/180,Ygrid(1)*pi/180,Xgrid(1)*pi/180));

%% Create Output Directory

if exist(output_Dir,'dir')==0
    mkdir(output_Dir)
end

%% Read Background Seismicity

if ~use_background
    mu_xy_Ml = zeros(length(Ycgrid),length(Xcgrid));
else
    mu_xy_Ml = load('background seismicity.txt');
end

%%% Background seismicity becomes a column, mu_xy_Ml_t is the transformed mu_xy_Ml
mu_xy_Ml_t = zeros(num_grid,1);
count=1;
for j=1:length(Ycgrid)
    for k=1:length(Xcgrid)
        mu_xy_Ml_t(count) = mu_xy_Ml(j,k);
        count=count+1;
    end
end

Nbo = sum(sum(mu_xy_Ml))*dA*T_start;   

%%  MCMC algorithm 

if do_MCMC_updating == 1   

    display(' ')
    display('---------- Use MCMC to find the parameters of ETAS model')   
    
    %%% Define the seq       
    indexSeq = find(time_MS >= 0 & time < T_start & M >= Mc); 
    
    %%% Calculate r
    r = cell(length(indexSeq)-1,1);
    for j=2:length(indexSeq)
        Grdata  = [latitude(indexSeq(j)),longitude(indexSeq(j))]*pi/180;           % Read the aftershock latitude and longitude, gradi a radian
        Grpdata = [latitude(indexSeq(1:j-1)),longitude(indexSeq(1:j-1))]*pi/180;   % Read the aftershock latitude and longitude, gradi a radian    
        xv = topgeo(Grpdata(:,1),Grpdata(:,2),Grdata(1),Grdata(2));            
        r{j-1,1} = sqrt(xv(:,1).^2+xv(:,2).^2);
    end
    
    %%% Find mu_xiyi_Ml
    mu_xiyi_Ml = zeros(1,length(indexSeq));
    for j=1:length(indexSeq)
        [~,indx_lat] = min(abs(Ycgrid-latitude(indexSeq(j))));
        [~,indx_lon] = min(abs(Xcgrid-longitude(indexSeq(j))));
        mu_xiyi_Ml(j) = mu_xy_Ml(indx_lat,indx_lon);
    end
        
    %%% Main Loop (for MCMC)
    samples = sample_posterior_MCMC_script_ETAS_spatialModel_new_Blockwise(indexSeq,M,time_MS,tstart,Mc,r,rxy,dA,do_calculate_K,consider_approximate,Nbo,mu_xiyi_Ml,priorPDF_par,THETA);
    
       
    %%% save
    filename = fullfile(output_Dir,'\samples.mat');
    save(filename,'samples','indexSeq','-v7.3');
    
else
        
    %%% load
    load(fullfile(output_Dir,'\samples.mat'));

end    
    
%% Finding the Robust Estimate 

if do_find_Robust_estimate == 1

display(' ')     
display('---------- Perform Robust estimation for the number of events') 
display(' ')     

%%% Assignments
sampleN_mgrM = zeros(num_grid,size(samples,2),length(vecM));

M_seqgen    = cell(1,size(samples,2));
time_seqgen = cell(1,size(samples,2));
Lon_seqgen  = cell(1,size(samples,2));
Lat_seqgen  = cell(1,size(samples,2));
r_seqgen    = cell(1,size(samples,2));
Int_t_lambdaseq  = cell(1,size(samples,2));

%%% Main Loop
for j=1:size(samples,2)
    
    display(['             Generated Seq no. ',num2str(j)]) 
    
    %%% Generate seqgen
    [M_seqgen{1,j},time_seqgen{1,j},Lon_seqgen{1,j},Lat_seqgen{1,j},r_seqgen{1,j},Int_t_lambdaseq{1,j}] = generateSEQ (M(indexSeq),time_MS(indexSeq),rxy(indexSeq,:),tstart,tend,Mc,samples(:,j),Xcgrid,Ycgrid,Ggrid,Mmax,dA,consider_approximate,mu_xy_Ml_t);
    
    %%% Calculate Robust N
    display('                  Calculate N ...')
    display(' ')
    for k=1:length(vecM)
        for ngrid=1:num_grid            
            if ~isempty(M_seqgen{1,j})
                sampleN_mgrM(ngrid,j,k) = calculate_N([M(indexSeq);M_seqgen{1,j}],[time_MS(indexSeq);time_seqgen{1,j}],[rxy(indexSeq,ngrid);r_seqgen{1,j}(:,ngrid)],vecM(k),tstart,tend,Mc,samples(:,j),Int_t_lambdaseq{1,j}(:,ngrid),dA,mu_xy_Ml_t(ngrid));
            else
                sampleN_mgrM(ngrid,j,k) = calculate_N(M(indexSeq),time_MS(indexSeq),rxy(indexSeq,ngrid),vecM(k),tstart,tend,Mc,samples(:,j),[],dA,mu_xy_Ml_t(ngrid));
            end
        end
    end
    
end

filename = fullfile(output_Dir,'\robust estimate N.mat');
save(filename,'sampleN_mgrM','-v7.3')

end

%% Post-Processing

if do_N_test_S_test == 1  

if ~do_find_Robust_estimate
    filename = fullfile(output_Dir,'\robust estimate N.mat');
    load(filename)
end

%%% Calculate the Total Number of Events
N_robust_mean = zeros(length(vecM),1);
N_robust_p50  = zeros(length(vecM),1);
N_robust_p16  = zeros(length(vecM),1);
N_robust_p84  = zeros(length(vecM),1);
N_robust_p02  = zeros(length(vecM),1);
N_robust_p98  = zeros(length(vecM),1);
N_observed    = zeros(length(vecM),1);

P_mean = zeros(length(vecM),1);
P_p98  = zeros(length(vecM),1);

N_robust_mean_xy = zeros(length(Ycgrid),length(Xcgrid),length(vecM));
N_robust_p50_xy  = zeros(length(Ycgrid),length(Xcgrid),length(vecM));
N_robust_p16_xy  = zeros(length(Ycgrid),length(Xcgrid),length(vecM));
N_robust_p84_xy  = zeros(length(Ycgrid),length(Xcgrid),length(vecM));
N_robust_p02_xy  = zeros(length(Ycgrid),length(Xcgrid),length(vecM));
N_robust_p98_xy  = zeros(length(Ycgrid),length(Xcgrid),length(vecM));

for k=1:length(vecM)
    
    display(' ')     
    display(['             Calculate the robust estimate for the number of event with M >= ',num2str(vecM(k))])  
    display(' ')     
    
    %%% Calculate N_robust
    [N_robust_mean(k),N_robust_p50(k),N_robust_p16(k),N_robust_p84(k),N_robust_p02(k),N_robust_p98(k)] = ordered_statistic(sum(sampleN_mgrM(:,:,k))); 
    [N_robust_mean_grid,N_robust_p50_grid,N_robust_p16_grid,N_robust_p84_grid,N_robust_p02_grid,N_robust_p98_grid] = ordered_statistic(sampleN_mgrM(:,:,k));
    
    if k==1
        Nforexy = N_robust_mean_grid;
    end       
    
    N_robust_mean_xy(:,:,k) = (reshape(N_robust_mean_grid,length(Xcgrid),length(Ycgrid)))';
    N_robust_p50_xy(:,:,k)  = (reshape(N_robust_p50_grid ,length(Xcgrid),length(Ycgrid)))';
    N_robust_p16_xy(:,:,k)  = (reshape(N_robust_p16_grid ,length(Xcgrid),length(Ycgrid)))';
    N_robust_p84_xy(:,:,k)  = (reshape(N_robust_p84_grid ,length(Xcgrid),length(Ycgrid)))';
    N_robust_p02_xy(:,:,k)  = (reshape(N_robust_p02_grid ,length(Xcgrid),length(Ycgrid)))';
    N_robust_p98_xy(:,:,k)  = (reshape(N_robust_p98_grid ,length(Xcgrid),length(Ycgrid)))';
    
    N_observed(k) = length(find(time >= T_start & time < T_end & M >= vecM(k)));
    
    %%% Calculate the probability
    P_mean(k) = 1-exp(-N_robust_mean(k));
    
end

%%% S-test
Nobsxy = zeros(length(Ggrid),1);
indexSeq_Tstart_tend = find(time >= T_start & time < T_end & M >= Mc);
for j=1:length(indexSeq_Tstart_tend)
    [~,indx_dummy] = min(rxy(indexSeq_Tstart_tend(j),:),[],2);
    Nobsxy(indx_dummy)=Nobsxy(indx_dummy)+1;
end

[Prob1_Stest,Prob_new_Stest] = S_test(Nobsxy,Nforexy,sampleN_mgrM(:,:,1));

end


%% END
        