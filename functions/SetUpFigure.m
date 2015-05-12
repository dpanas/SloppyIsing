%-----------------------------------------------------------------------------------%
% function: SetUpFigure
%           a short function to prepare figure handle and subplots,
%           optimized so that exported as pdf fits nicely into the thesis.
%           A few options to accommodate most common needs (e.g. single
%           long subplot for raster, or four square subplots for MEA
%           layout)
%           
%           WORK IN PROGRESS
%
% dependancy: --- (meant to work with SetPanelAppearance and SetPanelLabel)
%
% input:   - how many subplots (1,2, 3 or 4 atm);
%          - [optional] whether alternative arrangement to default is
%          wanted;
%
% output:  - figure handle;
%          - a cell of strings specifying the subplot position; each time
%          something is to be plotted, the string should be called through:
%          eval(subplot_cell{i})
%
%          DAP Mar 2015
%-----------------------------------------------------------------------------------%

function [fig_handle,subplot_cell]=SetUpFigure(how_many_subplots,arrangement)

fig_handle = figure;
subplot_cell = cell(how_many_subplots,1);
    
if how_many_subplots==1
    if ~exist('arrangement') 
        set(gcf,'Position',[230 200 920 430])
        subplot_cell{1} = 'subplot(''position'',[.27 .15 .5 .72])';
    else
        set(gcf,'Position',[230 200 920 350])
        subplot_cell{1} = 'subplot(''position'',[.09 .18 .88 .68])';
    end
elseif how_many_subplots==2
    if ~exist('arrangement') 
        set(gcf,'Position',[230 200 920 430])
        subplot_cell{1} = 'subplot(''position'',[.09 .15 .37 .72])';
        subplot_cell{2} = 'subplot(''position'',[.6 .15 .37 .72])';    
    else
        set(gcf,'Position',[230 100 920 700])
        subplot_cell{1} = 'subplot(''position'',[.09 .56 .88 .34])';
        subplot_cell{2} = 'subplot(''position'',[.09 .09 .88 .34])'; 
    end
elseif how_many_subplots==3
    if ~exist('arrangement') 
        set(gcf,'Position',[230 200 920 350])
        subplot_cell{1} = 'subplot(''position'',[.09 .18 .23 .68])';
        subplot_cell{2} = 'subplot(''position'',[.41 .18 .23 .68])';    
        subplot_cell{3} = 'subplot(''position'',[.73 .18 .23 .68])';
    else
        set(gcf,'Position',[230 150 920 640])
        subplot_cell{1} = 'subplot(''position'',[.09 .7 .88 .23])';
        subplot_cell{2} = 'subplot(''position'',[.09 .4 .88 .23])';    
        subplot_cell{3} = 'subplot(''position'',[.09 .1 .88 .23])';   
    end
elseif how_many_subplots==4
    if ~exist('arrangement') 
        set(gcf,'Position',[230 10 920 840])
        subplot_cell{1} = 'subplot(''position'',[.09 .58 .37 .36])';
        subplot_cell{2} = 'subplot(''position'',[.6 .58 .37 .36])';    
        subplot_cell{3} = 'subplot(''position'',[.09 .1 .37 .36])';   
        subplot_cell{4} = 'subplot(''position'',[.6 .1 .37 .36])';   
    else
        set(gcf,'Position',[230 200 920 300])
        subplot_cell{1} = 'subplot(''position'',[.04 .21 .17 .57])';
        subplot_cell{2} = 'subplot(''position'',[.27 .21 .17 .57])';    
        subplot_cell{3} = 'subplot(''position'',[.5 .21 .17 .57])';   
        subplot_cell{4} = 'subplot(''position'',[.73 .21 .17 .57])';         
    end
elseif how_many_subplots==5
    set(gcf,'Position',[230 10 920 840])
    subplot_cell{1} = 'subplot(''position'',[.09 .81 .88 .14])';
    subplot_cell{2} = 'subplot(''position'',[.09 .63 .88 .14])';    
    subplot_cell{3} = 'subplot(''position'',[.09 .45 .88 .14])';   
    subplot_cell{4} = 'subplot(''position'',[.09 .27 .88 .14])'; 
    subplot_cell{5} = 'subplot(''position'',[.09 .09 .88 .14])'; 
else
    disp('Too many subplots wanted')
end

end

