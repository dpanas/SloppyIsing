%-------------------------------%
% function: CalculateModelStats
%           calculates the first and second order statistics of the data,
%           under assumption of Ising model, given model parameters
%
% dependancy: ---
%
% input:   - 'fields' of Ising model for each channel;
%          - 'interactions' of Ising model between each pair of channels;
%
% output:  - mean spikes in each channel predicted by model;
%          - mean crosspikes between channels predict. by model;
%          - entropy of the prob distribution defined by the model;
%
% DAP April 2013
% !!! for now v basic, no user control of algorithm (no Monte Carlo) !!!
% !!! no error control !!! -1 and 1 coding assumed
%-------------------------------%


function [s_ising,ss_ising,E_ising] = CalculateModelStats(H,J)

N = size(J,2);
s_ising = zeros(1,N);
ss_ising = zeros(N,N);
Z = 0;
E_ising = 0; 

% exact algorithm: go through all possible state vectors and calculate the
% partition function and, using that, model statistics (first and second
% order)
for i = 1:2^N
    
    % state vector:
    s = (-1)*ones(1,N);
    for k = 1:N
        if (mod(floor(i/(2^(N-k))),2)==0) s(k) = 1; end
    end   
    % unnormalized probability of the state vec. under model parameters:
    prob = 0;
    for l=1:N
        prob = prob + H(l)*s(l);
        for ll=(l+1):N
            prob = prob + J(l,ll)*s(l)*s(ll);
        end
    end
    prob = exp(prob);  
    E_ising = E_ising - prob*log2(prob);
    Z = Z + prob;
    s_ising = s_ising + s*prob;
    ss_ising = ss_ising + s'*s*prob;
    
% %     HERE IS THE ALTERNATIVE IMPLEMENTATION:
%     
%     k = N - (1:1:N);
%     s = 1-2*mod(floor(i./(2.^k)),2);
%     
%     prob = H*s' + s*(triu(J-eye(N))*s');
%     prob = H*s' + s*J*s';
%     prob = exp(prob);  
%     E_ising = E_ising - prob*log2(prob);
%     Z = Z + prob;
%     s_ising = s_ising + s*prob;
%     ss_ising = ss_ising + s'*s*prob;
    
    
end
s_ising = s_ising/Z;
ss_ising = ss_ising/Z;
E_ising = E_ising/Z + log2(Z) ;

end