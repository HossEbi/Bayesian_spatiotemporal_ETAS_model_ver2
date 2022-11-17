% Calculate the kernel PDF
% Seismicity forecasting based on a Bayesian spatio?temporal ETAS model
% written by: Hossein Ebrahimian 
% Last update: 11/2022

function pdfValue = kernelPDF(theta,seeds,type)

weights = calculateWeights(seeds,type);

pdfValue = calculateKernel(theta,seeds,weights);

return
