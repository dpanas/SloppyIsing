%-------------------------------%
% function: GetSameChannels
%           Filter out channels that aren't active across all recordings.
%
% dependancy: ---
%
% input:  - cell holding channel keys (the 'chankey' output of GetSpikes)
%           of all the recordings of interest;
%
% output:  - a single array with numbers and coordinates of channels that
%            are active across all the recordings of interest;
%          - optionally also there is a possibility to output chankeys for
%            each recording containing only common channels;
%
% DAP April 2013
% !!! no error control !!!
%-------------------------------%

function [chankey,chankeyAllNew] = GetSameChannels(chankeyAll)

nof = max(size(chankeyAll));

prop = chankeyAll{1}(:,1);
% go over recordings:
for i = 2:nof
    to_discard = [];
    % go over channels from the first recording and check if they appear in
    % following recordings - if not, discard:
    for j = 1:length(prop)
        %disp('bjsb')
        prop2 = find(chankeyAll{i}(:,1)==prop(j));
        %length(prop2)
        if length(prop2)<1
            to_discard = [to_discard j];            
        end
    end
    prop(to_discard) = [];
end
prop2 = [];
for i = 1:length(prop)
    prop2 = [prop2 find(chankeyAll{1}(:,1)==prop(i))];
end
chankey = chankeyAll{1}(prop2,1:3);

chankeyAllNew = cell(nof,1);
for i = 1:nof
   prop2 = [];
   for j = 1:length(prop)
       prop2 = [prop2 find(chankeyAll{i}(:,1)==prop(j))];
   end
   chankeyAllNew{i} = chankeyAll{i}(prop2,:);
end

end