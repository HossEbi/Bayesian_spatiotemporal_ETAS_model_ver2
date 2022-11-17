% Function for generating SEQgen
% Seismicity forecasting based on a Bayesian spatio?temporal ETAS model
% written by: Hossein Ebrahimian 
% Last update: 11/2022

function [Mseq,timeseq,Lonseq,Latseq,rseq,Int_t_lambdaseq] = generateSEQ (Mi, Ti, ri, tstart, tend, Ml, sampletheta, Xgrid, Ygrid, Ggrid, Mmax, dA, tgen_approximate, mu_xy_Ml)

%% Initial Assignments

Mseq    = [];
timeseq = [];
Lonseq  = [];
Latseq  = []; 
rseq    = []; 
Int_t_lambdaseq = [];

mub = sum(mu_xy_Ml)*dA;

%% Generating 1st Magnitude and time

%%% Estimate the integral for distance
if ~tgen_approximate
    Iri =  calculate_Ir (ri, sampletheta(6:end), dA, Mi);
else
    Iri = [];
end

%%% generate the 1st set of data (magnitude and time)
[Mgen,tgen,numRepeat] = generate_m_time (Mi, Ti, Iri, tstart, tend, Ml, sampletheta, Mmax, mub);

%% Main loop for generating other parameters

if tgen > tend
    
    display(['                  tgen = ',num2str(tgen),' > T_end = ',num2str(tend),', NO sequence will be generated!!!'])
    
else
    
    %%% generate the 1st set of lat and lon
    [rgen,Longen,Latgen,Int_t_lambda] = generate_R (Mi, Ti, ri, Mgen, tstart, tgen, Ml, sampletheta, Xgrid, Ygrid, Ggrid, mu_xy_Ml);

    count = 0;
    while tgen <= tend
        
        display(['                  tgen = ',num2str(tgen),' / ',num2str(tend),' , Mgen = ',num2str(Mgen),' - Thinning iteration = ',num2str(numRepeat)])
   
        count = count + 1;
        Mseq(count,:)    = Mgen;
        timeseq(count,:) = tgen;
        Lonseq(count,:)  = Longen;
        Latseq(count,:)  = Latgen;
        rseq(count,:)    = rgen;
        Int_t_lambdaseq(count,:) = Int_t_lambda;
        
        if ~tgen_approximate
            Iri =  calculate_Ir ([ri;rseq], sampletheta(6:end), dA, [Mi;Mseq]);
        end
    
        [Mgen,tgen,numRepeat] = generate_m_time ([Mi;Mseq], [Ti;timeseq], Iri, tgen, tend, Ml, sampletheta, Mmax, mub);
        [rgen,Longen,Latgen,Int_t_lambda] = generate_R ([Mi;Mseq], [Ti;timeseq], [ri;rseq], Mgen, timeseq(end), tgen, Ml, sampletheta, Xgrid, Ygrid, Ggrid, mu_xy_Ml);
        
    end

end
 
end




    
    
    
    
    
    
    
    