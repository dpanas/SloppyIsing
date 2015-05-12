%-------------------------------%
% function: GetDataRates
%           Get the distribution of probabilities of each pattern of spikes
%           for provided spike matrix of binned spikes.
%
% dependancy: ---
%             !!! -1 and 1 coding assumed !!!
%
% input:    - matrix with binned channels (column = channel);
%
% output:  - matrix with data rates; (all per bin)
%
% !!! for now v basic, just rate of states by model vs data, exact count !!!
% !!! no error control !!!
%-------------------------------%

function [data_rates,data_zipf] = GetDataRates(data)

N = size(data,2);
bins = size(data,1);
data_rates = zeros(1,2^N);
pat_len = zeros(1,2^N);

% exact algorithm: go through all possible state vectors and calculate the
% probabilities -> and rates from them; also rates from the data
for i = 1:2^N
    % state vector:
    s = (-1)*ones(1,N);
    for k = 1:N
        if (mod(floor(i/(2^(N-k))),2)==0) s(k) = 1; end
    end
    pat_len(i) = sum((s+1)./2);
    data2 = data;
    for l=1:N
       data2(:,l) = data2(:,l) - s(l);
       data2 = data2(find(data2(:,l)==0),:);
    end
    data_rates(i) = size(data2,1);
end
data_rates = data_rates/bins;

data_zipf = zeros(1,N+1);
for j=1:N
    prop = find(pat_len==j);
    data_zipf(j+1) = sum(data_rates(prop));
end
prop = find(pat_len==0);
data_zipf(1)=sum(data_rates(prop));

end