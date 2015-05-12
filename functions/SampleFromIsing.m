%-------------------------------%
% function: SampleFromIsing
%           a discretized version of the Ising model distribution of
%           firing patterns (as if data was according to Ising, but due to
%           limited samples the prob distribution is discretized, some
%           patterns don't appear etc) + plus if needed a surrogate dataset
%
% dependancy: ---
%
% input:  - Ising model rate of pattern occurence;
%         - number of desired time bins of the data;
%         - additional variable indicating whether to generate data or not;
%
% output:  - sampled rate of pattern occurence;  
%          - surrogate spiking data (already binned and -1 and 1 trains);             
%
% !!! no error control !!!
%-------------------------------%

function [test_rate,test_spikes]=SampleFromIsing(model_rate,Nbins,don)

test_rate = zeros(1,length(model_rate));
% the first pattern prob. calculated separately:
test = rand(Nbins,1);   % a random vector of the data length
prop = find(test<=model_rate(1));   % counting how many random numbers fall within a desired range
test_rate(1) = length(prop);
% and all the following, similarly:
for j=2:(length(model_rate))        
    prop = test(find(test>sum(model_rate(1:j-1))));
    prop = prop(find(prop<=sum(model_rate(1:j))));
    test_rate(j) = length(prop);
end

N = log2(length(test_rate));
test_spikes = [];
% if more than two inputs, then surrogate data is generated:
if nargin>2
    for i = 1:2^N
        % state vector:
        s = (-1)*ones(1,N);
        for k = 1:N
            if (mod(floor(i/(2^(N-k))),2)==0) s(k) = 1; end
        end
        test_spikes = [test_spikes; meshgrid(s,ones(1,test_rate(i)))];
    end
    % a random permutation of bin numbers to assign patterns to:
    test_spikes = test_spikes(randperm(Nbins),:);
end

test_rate = test_rate./Nbins;  % normalizing the counts
end