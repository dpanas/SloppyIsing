%-------------------------------%
% function: FitModelIndep
%           Fits an independent model to the data - for comparison with
%           Ising model; it's simply a calculation of possible pattern from
%           the rates.
%
% dependancy: -
%             !!! -1 and 1 coding assumed !!!
%
% input:  - mean spikes in each channel;
%
% output:  - 'mean fields' of independent model for each channel;
%          - entropy of the independent model distribution;
%
% DAP August 2013
% !!! no error control !!! 
%-------------------------------%


function [rate_indep,E_indep] = FitModelIndep(s_data)

N = size(s_data,2);
rate_indep = ones(1,2^N);
E_indep = 0;

% converting mean spike marginals into probabilities of firing:
s_data = (s_data+1)./2;

for i = 1:2^N
    % state vector:
    s = (-1)*ones(1,N);
    for k = 1:N
        if (mod(floor(i/(2^(N-k))),2)==0) s(k) = 1; end
    end
    % probability of the state vector:
    prob = s.*s_data - ((s-1)./2);
    for l=1:N
        rate_indep(i) = rate_indep(i)*prob(l);
    end
    E_indep = E_indep - rate_indep(i)*log2(rate_indep(i));
end

end