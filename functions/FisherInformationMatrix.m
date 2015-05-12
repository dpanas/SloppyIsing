%-------------------------------%
% function: FisherInformationMatrix
%           calculates the Fisher Information Matrix for a single Ising
%           model, see Machta and Sethna publications
%
% dependancy: ---
%
% input:   - probability distribution of the model;
%          - first order marginals according to the model;
%          - second order marginals according to the model;
%
% output:  - 
%
% DAP April 2014
% !!! for now v basic, no user control of algorithm (no Monte Carlo) !!!
% !!! no error control !!! -1 and 1 coding assumed
%-------------------------------%

function [FIM]=FisherInformationMatrix(model_rate,s_ising,ss_ising)

N = size(ss_ising,2);
FIM = zeros(N*(N+1)/2,N*(N+1)/2);
param = s_ising;
for j=1:N
    param = [param ss_ising(j,j+1:end)];
end

for i = 1:2^N
    
    % state vector:
    s = (-1)*ones(1,N);
    for k = 1:N
        if (mod(floor(i/(2^(N-k))),2)==0) s(k) = 1; end
    end   
    SS = s'*s;
    for j=1:N
        s = [s SS(j,j+1:end)];
    end
    
    FIM = FIM + model_rate(i).*(s-param)'*(s-param);
    
end

end