%-------------------------------%
% function: GetCorrelations
%           Calculates the Pearson correlation coefficient between each
%           pair of channels listed in 'chankey' (channels treated as binary
%           vectors of spiking / not spiking - requires the spiketrains to
%           be binned already). 
%           Also provides a statistical test for finding common spikes
%           between each pair of channels - can be used later to assess
%           significance of correlations (recommended with multiple
%           comparisons correction).
%           Also calculates distance between each pair of channels.
%
% dependancy: ---
%
% input:   - channel key: [active channel ID; x coordinate on MEA; y
%           coordinate on MEA; (...)]; 
%          - a cell array; each cell is a vector with indices of bins that 
%           hold spikes; (spikes need to be binned previously, ID 
%           corresponds to first column of chankey)
%          - number of bins;
%
% output:  - matrix of correlations: (i,j) between i'th and j'th channel as
%            indexed in the 'chankey';
%          - matrix of distances: (i,j) between i'th and j'th channel as
%            indexed in the 'chankey';
%          - matrix of CDF of hypergeometric distribution with parameters 
%            corresponding to i'th and j'th channel: basically the
%            probability of randomly finding the observed number of common spikes
%            or fewer, given firing rates and number of samples (bins)
%
% DAP Mar 2014
% !!! no error control !!!
%-------------------------------%

function [correlations,distances,probtest]=GetCorrelations(chankey_binned,spikes_binned,Nbins)

dims = size(chankey_binned,1);
correlations = zeros(dims,dims);
distances = zeros(dims,dims);
probtest = zeros(dims,dims);

for i=1:dims;
    a=-1*ones(1,Nbins);  
    a(spikes_binned{chankey_binned(i,1)})=1;
    Na = sum((a+1)./2);
    ma = 2*Na./Nbins - 1;
    sa = std(a);
    
    for j=i+1:dims;
        b=-1*ones(1,Nbins);             
        b(spikes_binned{chankey_binned(j,1)})=1;
        Nb = sum((b+1)./2);
        mb = 2*Nb./Nbins - 1;
        sb = std(b);
        C = ((a+1)./2)*((b+1)./2)';
        clear b
        probtest(i,j) = hygecdf(C,Nbins,Na,Nb);
        correlations(i,j) = (Nbins+4*C-2*Na-2*Nb-Nbins*ma*mb)./sa/sb/(Nbins-1);
        distances(i,j) = sqrt( ( chankey_binned(i,2)-chankey_binned(j,2) )^2 + ( chankey_binned(i,3)-chankey_binned(j,3) )^2 );        
    end
    
end


end