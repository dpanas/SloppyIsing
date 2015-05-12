%-------------------------------%
% function: InformationMeasures
%           calculates entropies of given distributions and Shannon-Jensen
%           divergence between the distributions
%
% dependancy: ---
%
% input:   - probability distribution 1;
%          - probability distribution 2; (assumed that these are over same
%          space of states obviously)
%
% output:  - entropy of distribution 1 [bits per bin];
%          - entropy of distribution 2 [bits per bin];
%          - Kullback-Leibler divergence [bits];
%          - Shannon-Jensen divergence [bits];
%
% !!! no error control !!!
%-------------------------------%

function [S1,S2,KL,SJ] = InformationMeasures(distr1,distr2)


dd = (distr1 + distr2) ./2;
prop1 = find(distr1>0);
prop2 = find(distr2>0);
SJ = ( sum(distr1(prop1).*log2(distr1(prop1)./dd(prop1))) )./2 + ( sum(distr2(prop2).*log2(distr2(prop2)./dd(prop2))) )./2;
KL = sum(distr1(prop1).*log2(distr1(prop1)./distr2(prop1)));

distr1 = distr1(prop1);
distr2 = distr2(prop2);
S1 = -sum(distr1.*log2(distr1));
S2 = -sum(distr2.*log2(distr2));

end