% Function for calculating the spatial distribution of the generated events
% Seismicity forecasting based on a Bayesian spatio?temporal ETAS model
% written by: Hossein Ebrahimian 
% Last update: 11/2022

function [rgen,Longen,Latgen,Int_t_lambda] = generate_R(Mi, Ti, ri, Mgen, tstart, tgen, Ml, theta, Xgrid, Ygrid, Ggrid, mu_xy_Ml)

pxy = zeros(length(Ggrid),1);

Int_t_lambda = zeros(1,length(Ggrid));

for j=1:length(Ggrid)
    [pxy(j),Int_t_lambda(j)] = calculate_pxy(Mi, Ti, ri(:,j), Mgen, tstart, tgen, Ml, theta, mu_xy_Ml(j));
end

pxy = (reshape(pxy,length(Xgrid),length(Ygrid)))'/sum(pxy);

[Longen,Latgen] = sample_pxy(pxy,Xgrid,Ygrid);

rgen = calculate_rxy(Latgen,Longen,Ggrid);

end