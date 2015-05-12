%-------------------------------%
% function: GetBinnedStats
%           Calculates stats in binned data. For the moment only absolute
%           mean and variance.
%
% dependancy: ---
%
% input:   - variable Y that is to provide the stats;
%          - variable X that defines the bins;
%          - bin edges for the binning variable X;
%
% output:  - bin centers (for plotting, if edges are not spaced evenly the
%            bin centers won't be either!);
%          - the mean absolute value of variable Y in each bin of X;
%          - the variance of variable Y in each bin of X;
%
% DAP Feb 2015
% !!! no error control !!!
%-------------------------------%


function [bin_centers,bin_mean,bin_variance] = GetBinnedStats(variable,variable_for_bins,bin_edges)

% set any NaNs to zero:
prop = find(isnan(variable));
variable(prop) = 0;
% bin centers for plotting:
bin_centers = bin_edges(1:end-1) + diff(bin_edges)./2;
bin_mean = [];
bin_variance = [];
[~,bin] = histc(variable_for_bins,bin_edges);
% ignoring the last bin, as the last counts only values ==bin_edges(end)
for i=1:length(bin_edges)-1;
    prop = find(bin==i);
    bin_mean = [bin_mean  mean(abs(variable(prop)))];
    bin_variance = [bin_variance var(variable(prop))];
end

end