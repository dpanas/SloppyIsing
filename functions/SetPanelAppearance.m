%-------------------------------%
% function: SetPanelAppearance
%           a short function to prepare plots for being of 'publishable'
%           quality - fixing the title, axes ranges, labels; pretty much
%           all inputs are optional, i.e. can be left [] to leave blank
%
% dependancy: --- 
%
% input:   - string with the title;
%          - string for x axis label;
%          - range for the x axis;
%          - range for the x axis grid, if grid is to be displayed;
%          - string for x axis label;
%          - range for the x axis;
%          - range for the x axis grid, if grid is to be displayed;
%
% DAP Dec 2014
%-------------------------------%

function SetPanelAppearance(string_title,string_xaxis,range_xaxis,grid_xaxis,string_yaxis,range_yaxis,grid_yaxis)

% if needed, title of the plot:
if ~isempty(string_title)
    title_handle = title(string_title);
    set(title_handle,...
        'FontName','Arial',...
        'FontSize',19,...
        'color',[.25 .25 .25])    
end

% set up the X axis:
if ~isempty(string_xaxis)
    xlabel_handle = xlabel(string_xaxis);
    set(xlabel_handle,...
        'FontName','Arial',...
        'FontSize',19,... 
        'color',[.25 .25 .25])
end
if ~isempty(range_xaxis)
    xlim([min(range_xaxis) max(range_xaxis)])
    if length(range_xaxis>2)
        set(gca,'XTick',range_xaxis)
    end
end

% set up the Y axis:
if ~isempty(string_yaxis)
    ylabel_handle = ylabel(string_yaxis);
    set(ylabel_handle,...
        'FontName','Arial',...
        'FontSize',19,... 
        'color',[.25 .25 .25])
end
if ~isempty(range_yaxis)
    ylim([min(range_yaxis) max(range_yaxis)])
    if length(range_yaxis)>2
        set(gca,'YTick',range_yaxis)
    end
end

% set up general properties of plot and axes:
set(gca,...
    'Box','off',...
    'color',[1 1 1],...
    'Xcolor',[.25 .25 .25],'Ycolor',[.25 .25 .25],...
    'LineWidth',2,...
    'TickDir','in',...
    'TickLength',[.02 .02],...
    'FontName','Arial',...
    'FontSize',16)

% set up gridlines - manually, so that they look better (and are not
% necessarily at the same locations as ticks):
if ~isempty(grid_xaxis)
    if length(grid_xaxis)<2
        grid_xaxis=range_xaxis;
    end
    if length(grid_xaxis)<2
        disp('Cannot set up grid with just one grid line')
    else
        ygrid = repmat([min(range_yaxis); max(range_yaxis)],1,numel(grid_xaxis));
        xgrid = [grid_xaxis; grid_xaxis];
        hold on
        plot(xgrid,ygrid,'color',[.5 .5 .5],'lineWidth',.5)
    end
end
if ~isempty(grid_yaxis)
    if length(grid_yaxis)<2
        grid_yaxis=range_yaxis;
    end
    if length(grid_yaxis)<2
        disp('Cannot set up grid with just one grid line')
    else
        xgrid = repmat([min(range_xaxis); max(range_xaxis)],1,numel(grid_yaxis));
        ygrid = [grid_yaxis; grid_yaxis];
        hold on
        plot(xgrid,ygrid,'color',[.5 .5 .5],'lineWidth',.5)
    end
end

end
