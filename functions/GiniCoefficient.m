%-------------------------------%
% function: GiniCoefficient
%           calculates the Gini coefficient for a data vector - that is,
%           measure of sparsity of the distribution of values in the
%           vector; absolute values are taken into account!
%
% dependancy: ---
%
% input:   - vector of data points;
%
% output:   - value of the Gini coefficient
%
% DAP July 2014
%-------------------------------%


function [gini] = GiniCoefficient(data)

data = sort(abs(data));
data = cumsum(data);
data = data./(max(data));
gini = 1 - 2*(sum(data./length(data))-1./length(data));

end