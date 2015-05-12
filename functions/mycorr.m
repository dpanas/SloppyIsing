%-------------------------------%
% function: mycorr
%           just a simple function for the expected value of pairwise
%           interactions, which is a stats used for computing Ising model
%           fits. Also, alternatively (when additional input is provided)
%           computes Pearson correlation coefficient (to be independent of
%           statistical toolbox if needed).
%
% dependancy: ---
%
% input:  - data matrix (each vector is a column, dat(:,i));
%         - additional input if the result is supposed to be pearson 
%           correlation coefficient and not pairwise interaction;
%
% DAP April 2013
% !!! no error control !!!
%-------------------------------%

function c = mycorr(dat,pearson)
[m,n]=size(dat);
c = zeros(n,n);
for i=1:n
    for j=1:n
        c(i,j) = sum((dat(:,i)).*(dat(:,j)))/(m);
        if nargin>1
            c(i,j) = sum((dat(:,i)-mean(dat(:,i))).*(dat(:,j)-mean(dat(:,j))))/std(dat(:,i))/std(dat(:,j))/(m-1);
        end
    end
end

end