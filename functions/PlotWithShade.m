%-------------------------------%
% function: PlotWithShade
%           a short function to plot semi-transparent shaded areas, e.g. to
%           use for showing confidence intervals or std around a mean.
%
% dependancy: --- 
%
% input:   - x vector; (can be empty one [])
%          - main line to be plotted; (can be an empty one [])
%          - line of the lower boundary for shaded region;
%          - line of the upper boundary for shaded region;
%          - color vector (rgb); (can be empty one [] or skipped altogether with transparency)
%          - transparency (on the scale from 0 to 1); (can be skipped)
% output:  - returns handle to the filled regions between two lines;
%          - returns handle to the overlaid main curve, if it is provided;
%
% DAP Dec 2014
%-------------------------------%

function [handle1,handle2] = PlotWithShade(xline,yline,ylineminus,ylineplus,col,transparency);

X1 = 1:1:length(ylineminus);
X2 = length(ylineminus):-1:1;

% setting a default color if not supplied:
if ~exist('col')
    col = [.8 .2 .4];
elseif isempty(col)
    col = [.8 .2 .4];
elseif (length(col)~=3)
    disp('Input argument for color should be a vector length 3')
end

% setting a default transparency if not supplied:
if ~exist('transparency')
    transparency = .4;
elseif isempty(transparency)
    transparency = .4;
end

if ~isempty(xline)
    handle1 = fill([[xline(X1)] [xline(X2)]],[[ylineminus(X1)] [ylineplus(X2)]],col,'FaceAlpha',transparency);
    set(handle1,'EdgeColor','None')
    if ~isempty(yline)
        hold on
        handle2 = plot(xline,yline,'lineWidth',2,'color',col);
    end
else
    handle1 = fill([[X1] [X2]],[[ylineminus(X1)] [ylineplus(X2)]],col,'FaceAlpha',transparency);
    set(handle1,'EdgeColor','None')
    if ~isempty(yline)
        hold on
        handle2 = plot(yline,'lineWidth',2,'color',col);
    end
end

end