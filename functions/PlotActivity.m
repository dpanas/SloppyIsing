%-------------------------------%
% function: PlotActivity
%           plots mean log activity, histogram of rates, a sorted raster 
%           and histogram of spike counts in 5ms time bins, and
%           a surface plot of mean activity
%
% dependancy: ---
%
% input:   - recording time in sec;
%          - channel key: [active channel no; x coordinate on MEA; y
%          coordinate on MEA];
%          - a cell array; each cell is a vector with spike times;

% !!! no error control !!!
%-------------------------------%

function PlotActivity(rtime,chankey,spikes)

% plot log mean activity:
figure
Wlen = 1400;
Whei = 450;
set(gcf,'Position',[100,300,Wlen,Whei])
subplot(121)
MeanAct0 = NaN(65,65);
for i=1:size(chankey,1)
    MeanAct0(chankey(i,2),chankey(i,3)) = log(chankey(i,4)/rtime);
end
cscale = [min(min(MeanAct0)) max(max(MeanAct0))];
subplot('position',[0.05 0.15 0.7*Whei./Wlen 0.7])
pcolor(MeanAct0')
shading flat
axis off
colormap jet
caxis(cscale)
title(['Total active channels: ',num2str(size(chankey,1))],'fontsize',13)
subplot('position',[(0.05+ 0.7*Whei./Wlen +0.1) 0.15 (1-(0.05+ 0.7*Whei./Wlen +0.1)-0.05 -0.05 )/2 0.7 ])
edges=0:0.002:max(chankey(:,4)./rtime);
H = histc(chankey(:,4)./rtime,edges);
bar(edges,H)
xlim([-0.1 max(chankey(:,4)./rtime)+0.1])
title(['Mean rate: ',num2str(mean(chankey(:,4)./rtime))],'fontsize',13)
xlabel('Rate [Hz]','fontsize',13)
ylabel('Counts','fontsize',13)
subplot('position',[(0.05+ 0.7*Whei./Wlen +0.1 +(1-(0.05+ 0.7*Whei./Wlen +0.1)-0.05 -0.05 )/2+ 0.05) 0.15 (1-(0.05+ 0.7*Whei./Wlen +0.1)-0.05 -0.05 )/2 0.7 ])
loglog(edges,H,'r.')
xlabel('Log rate [Hz]','fontsize',13)
ylabel('Log counts','fontsize',13)

% sort the channels according to mean rate:
chky = chankey;
chky = sortrows(chky,[4]);

spike_count = [];              % spike count for the histogram
isis = [];                     % interspike intervals from each channel pooled together
MeanAct = NaN(64*3,64*3);      % matrix for mean activation
vec     = zeros(64*3,1);       % vector to use for plotting MeanAct

figure
set(gcf,'Position',[50,200,1500,600])
subplot(211)
hold on
% go through sorted channels and plot raster, and record MeanAct
for i=1:length(chky)
    plot(spikes{chky(i,1)},i*ones(length(spikes{chky(i,1)}),1),'k.','MarkerSize',1)
    spike_count = [spike_count spikes{chky(i,1)}'];
    isis = [isis diff(spikes{chky(i,1)})'];
    x = chky(i,2);
    y = chky(i,3);
    MeanAct(((x-1)*3+1):x*3,((y-1)*3+1):y*3) = (chky(i,4)/rtime);
end
xlabel('Time [sec]','fontsize',13)
ylabel('Channel nr (sorted)','fontsize',13)
% plot the histrogram of spikes in 5ms time bins:
subplot(212)
hbin = 0.005;
edges=0:hbin:rtime;
H = histc(spike_count,edges);
bar(edges+hbin./2,H)
% info scales etc
xlim([0 rtime])
xlabel('Time [sec]','fontsize',13)
ylabel('Counts','fontsize',13)

figure
set(gcf,'Position',[200,100,1080,620])
edges = 0:0.0005:1;
h = histc(isis,edges);
bar(edges,h)
xlim([-0.05 0.95])
xlabel('Interspike interval [sec]','fontsize',13)
ylabel('Counts','fontsize',13)
title('ISI distribution','fontsize',13)
% figure
% for x=1:64
%     vec(((x-1)*3+1):x*3) = [x-0.49 x x+0.49];
% end
% surf(vec,vec,(MeanAct'))   %the x and y axis are swapped by surf, hence the comma
% % info scales etc
% title('Rates [Hz]','fontsize',13)
% xlim([0 65])
% xlabel('X coord','fontsize',13)
% ylim([0 65])
% ylabel('Y coord','fontsize',13)
    
end