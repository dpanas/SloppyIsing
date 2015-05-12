%-------------------------------%
% function: SetPanelLabel
%           a short function to add a label to the up and left of the panel
%
% dependancy: --- 
%
% input:   - label string;
%          - offset to left (fraction of panel size);
%          - offset up (fraction of panel size);
%
% DAP Mar 2015
%-------------------------------%

function SetPanelLabel(label,left,up)

aa = axis;
pos(1) = aa(1) - (aa(2)-aa(1))*left;
pos(2) = aa(4) + (aa(4)-aa(3))*up;
text(pos(1),pos(2),label,'fontName','Arial','fontsize',22,'fontweight','bold','color',[.15 .15 .15])

end