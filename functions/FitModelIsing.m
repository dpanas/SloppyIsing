%-------------------------------%
% function: FitModelIsing
%           Fits the Ising model to the data with the iterative scaling
%           algorithm (as in Tang et al. 2008); this is an exact algorithm
%           so not too many neurons are advised (less than 20)
%
% dependancy: - CalculateModelStats; - ApproximateParameters;
%             !!! -1 and 1 coding assumed !!!
%
% input:  - mean spikes in each channel;
%         - mean crosspikes between channels;
%         - learning rate (value between 0 and 1; 0.75 is fairly safe);
%         - max number of iterations;
%
% output:  - 'mean fields' of Ising model for each channel;
%          - 'interactions' of Ising model for each pair of channels;
%          - entropy of the Ising model distribution;
%          - number of iterations taken to converge;
% 
% DAP August 2013
% !!! no error control !!!
%-------------------------------%


function [H,J,e,c] = FitModelIsing(s_data,ss_data,miu,max_iter)


c = 0;
% initialize model parameters:
H = s_data;         
J = ss_data;
%[H,J] = ApproximateParameters(s_data,ss_data);
[s_model,ss_model,~]=CalculateModelStats(H,J);
% fitting the model until desired convergence:
while ((max(abs((s_data-s_model)./s_data))>0.00001) || (max(max(abs((ss_data-ss_model)./ss_data)))>0.00001) ) && (c<max_iter)
    H = H + miu*sign(s_data).*log(abs(s_data./s_model));
    J = J + miu*sign(ss_data).*log(abs(ss_data./ss_model));
    [s_model,ss_model,e]=CalculateModelStats(H,J);
    c = c+1;
end

end