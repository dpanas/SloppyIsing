%-------------------------------%
% function: R2Coefficient
%           calculates the coefficient of determination (r^2) between the
%           two provided vectors (which must be of same lengths and the
%           observations need to be paired / corresponding); the second is
%           treated as a function of the first
%
% dependancy: ---
%
% input:   - vector 1 of data points;
%          - vector 2 of data points;
%
% output:   - value of the coefficient of determination
%
% DAP July 2014
%-------------------------------%


function [r2] = R2Coefficient(data1,data2)

p = polyfit(data1,data2,1);
data_p = p(1)*data1 + p(2);
ss_tot = sum((data2-mean(data2)).^2);
ss_res = sum((data2-data_p).^2);
r2 = 1 - ss_res/ss_tot;

end