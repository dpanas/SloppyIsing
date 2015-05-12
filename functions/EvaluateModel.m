%-------------------------------%
% function: EvaluateModel
%           some measures to see how well model fits the data: atm only
%           comparison of the rates of state vectors from the data versus
%           those predicted by the model
%
% dependancy: ---
%             !!! -1 and 1 coding assumed !!!
%
% input:   - 'mean fields' of Ising model for each channel;
%          - interactions of Ising model;
%          - 'mean fields' of independent model for each channel
%          - matrix with binned channels (column = channel);
%
% output:  - matrix with data prediction rates; (all per bin)
%          - matrix with model prediction rates;
%          - matrix with prediction rates of independent (poisson) model;
%
% !!! for now v basic, just rate of states by model vs data, exact count !!!
% !!! no error control !!!
%-------------------------------%

function [data_rates,model_rates,indep_rates] = EvaluateModel(H,J,h,data)

N = size(J,2);
data_rates = zeros(1,2^N);
model_rates = zeros(1,2^N);
indep_rates = ones(1,2^N);

if nargin==4
    bins = size(data,1);
    % exact algorithm: go through all possible state vectors and calculate the
    % probabilities -> and rates from them; also rates from the data
    for i = 1:2^N
        % state vector:
        s = (-1)*ones(1,N);
        for k = 1:N
            if (mod(floor(i/(2^(N-k))),2)==0) s(k) = 1; end
        end
        % unnormalized probability of the state vec. under model parameters:
        prob = 0;
        prob2 = 0;
        data2 = data;
        for l=1:N
            prob = prob + H(l)*s(l);
            prob2 = prob2 + h(l)*s(l);
            for ll=(l+1):N
                prob = prob + J(l,ll)*s(l)*s(ll);
            end
           data2(:,l) = data2(:,l) - s(l);
           data2 = data2(find(data2(:,l)==0),:);
        end
        model_rates(i) = exp(prob);
        indep_rates(i) = exp(prob2);
        data_rates(i) = size(data2,1);
    end
    model_rates = model_rates/sum(model_rates);
    data_rates = data_rates/bins;
    indep_rates = indep_rates/sum(indep_rates);    
elseif nargin==2
    % exact algorithm: go through all possible state vectors and calculate the
    % probabilities -> and rates from them; also rates from the data
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
        model_rates(i) = exp(prob);
    end
    model_rates = model_rates/sum(model_rates);
end

end