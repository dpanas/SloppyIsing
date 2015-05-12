%-----------------------------------------------------------------------------------%
% function: PlotBoxplot
%           plots a boxplot of data provided in a cell (with optional
%           inputs of labels and colors for boxes); automatically detects
%           if data is of different lengths; style is fixed
%
% dependancy: ---
%
% input:   - a cell with arrays of values to be boxplotted;
%          - [optional] labels for each of the boxes (single cell of strings);
%          - [optional] colors for each of the boxes (cell of RGBs);
%
%           DAP Feb 2015
%-----------------------------------------------------------------------------------%

function handle_boxes = PlotBoxplot(values_cell,labels,colors)

values_boxplot = [];
grouping_var = [];
ss = size(values_cell);
for i=1:length(values_cell)
    if ss(1)>ss(2)
        values_cell{i} = values_cell{i}(:)';
    end
    values_boxplot = [values_boxplot values_cell{i}];
    grouping_var = [grouping_var i*ones(1,length(values_cell{i}))];
end

% re-plot boxplot:
if ~isempty(labels)
    bh = boxplot(values_boxplot,grouping_var,'width',0.66,'labels',labels,'labelorientation','inline');
else
    bh = boxplot(values_boxplot,grouping_var,'width',0.66);
end
% and get handle to provide for legend:
handle_boxes =  findobj(gca,'tag','Box');
h = findobj(gca,'Tag','Outliers');

% set the color and style of median mark and outliers:
for i=1:length(values_cell)
    set(bh(3,i),'color','w')
    set(bh(4,i),'color','w')
    set(bh(6,i),'lineWidth',2,'color','k')
    set(bh(4,i),'color','w')
    set(h(i),'Marker','+','MarkerSize',2,'MarkerEdgeColor',[.3 .3 .3])
end


if ~isempty(colors)
    % set the color of boxes and whisker:
    for i=1:length(values_cell)
        patch(get(handle_boxes(length(values_cell)+1-i),'XData'),get(handle_boxes(length(values_cell)+1-i),'YData'),'y','FaceColor',colors{i},'EdgeColor',colors{i},'FaceAlpha',.9);
        set(bh(5,i),'lineWidth',1,'color',colors{i})
        set(bh(1,i),'lineWidth',1,'linestyle','-','color',colors{i})
        set(bh(2,i),'lineWidth',1,'linestyle','-','color',colors{i})
    end
else
    % or alternatively just the width
    for i=1:length(values_cell)
        set(bh(5,i),'lineWidth',2,'color',[31./255 119./255 180./255])
        set(bh(1,i),'lineWidth',2,'linestyle','-','color','k')
        set(bh(2,i),'lineWidth',2,'linestyle','-','color','k')
    end
end
% and set font for the labels:
set(findobj(get(handle_boxes(1), 'parent'), 'type', 'text'),...
    'fontsize', 14,...
    'color', [.25 .25 .25] ,...
    'FontName','Arial')

end