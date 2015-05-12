% %-------------------------------%
% script:   Compare_RatesCorrelations
%           Script for some basic plots about the distributions of rates
%           and correlations, as well as comparison of those and
%           statistical tests between recordings.
%
% !!! This is a script, not a function - parameters need to be changed
% manually by the user upon each call !!!
%
% dependancy: PlotBoxplot, GetSameChannels
%
% DAP Apr 2015
% %-------------------------------%

% first clear the workspace:
clear

% variables are here:
%-------------------%
% names of files with recorded, filtered, and binned spikes:
% (same as for script_SampleRandomGroups and others)
filenames{1} = './results_spikes/chip136_0_spikes_new_prob10_bin5_filt';
filenames{2} = './results_spikes/chip136_1_spikes_new_prob10_bin5_filt';
filenames{3} = './results_spikes/chip136_2_spikes_new_prob10_bin5_filt';
filenames{4} = './results_spikes/chip136_4_spikes_new_prob10_bin5_filt';
filenames{5} = './results_spikes/chip136_5_spikes_new_prob10_bin5_filt';
filenames{6} = './results_spikes/chip136_6_spikes_new_prob10_bin5_filt';
%-------------------%

% (below no more variables, just code)
%-------------------%
close all
addpath('./functions/')
NOF = length(filenames);
chankeyAll=cell(NOF,1);
ratesAll=cell(NOF,1);
correlationsAll=cell(NOF,1);
distancesAll=cell(NOF,1);
rates_pool = [];
group_rates = [];
correlations_pool = [];
group_correlations = [];
% read in (and pool together) the rates and correlations of each recording:
cols = jet(NOF); % set up colors
disp(' ')
disp(['Script comparing the firing rates and correlations between recordings.'])
disp(' ')
for i=1:NOF    
    eval(['load ',filenames{i}])
    disp(['Loaded spike file: ',' '' ',filename_spikes,' '' '])
    chankeyAll{i} = chankey_binned;
    ratesAll{i} = chankey_binned(:,4)'./rtime;
    rates_pool = [rates_pool ratesAll{i}];
    group_rates = [group_rates i*ones(1,length(ratesAll{i}))];
    distancesAll{i} = distances(:);
    correlationsAll{i} = correlations(:);
    correlationsAll{i}=correlationsAll{i}(distancesAll{i}~=0);
    distancesAll{i}=distancesAll{i}(distancesAll{i}~=0);
    correlations_pool = [correlations_pool correlationsAll{i}'];
    group_correlations = [group_correlations i*ones(1,length(correlationsAll{i}))];
    col{i} = cols(i,:);
end
clear chankey chankey_binned spikes spikes_binned time_compute_corr i Info filename_rec filename_data correlations distances cols
disp(' ')
disp(' ')
% plots and stats of the distributions of rates:
figure(1)
set(gcf,'Position',[100,100,1080,720])
subplot('position',[0.07 0.57 0.38 0.33])
PlotBoxplot(ratesAll,[],col);
xlabel('Recording','fontsize',13)
ylabel('Firing rate [Hz]','fontsize',13)
title('BOXPLOT','fontsize',14)
subplot('position',[0.56 0.57 0.4 0.33])
hold on
string=['legend('];
p_kstest_all_rates = [];
for i=1:NOF
    prop=sort(ratesAll{i});
    plot(prop,[1:length(prop)]./length(prop),'color',col{i})
    string=[string,'''Rec ',num2str(i),''','];
    if i>1
        [~,p] = kstest2(ratesAll{i-1},ratesAll{i});
        p_kstest_all_rates = [p_kstest_all_rates p];
    end
end
string=string(1:end-1);
string=[string,',''location'',''Southeast'')'];
eval(string)
xlabel('Firing rate [Hz]','fontsize',13)
ylabel('Fraction of channels','fontsize',13)
title('CDF','fontsize',14)
clear prop string i
[p_kruskal_rates,anovatab_rates] = kruskalwallis(rates_pool,group_rates,'off');
disp(['Nonparametric ANOVA for the firing rates: ',num2str(p_kruskal_rates)])
disp('     and Kolmogorov-Smirnov test between consecutive recordings:  ')
vpa(p_kstest_all_rates,4)
disp(' ')
disp(' ')
% plots and stats for the distributions of correlations:
figure(2)
set(gcf,'Position',[100,100,1080,720])
subplot('position',[0.07 0.57 0.38 0.33])
PlotBoxplot(correlationsAll,[],col);
xlabel('Recording','fontsize',13)
ylabel('Correlation coefficient','fontsize',13)
title('BOXPLOT','fontsize',14)
subplot('position',[0.56 0.57 0.4 0.33])
hold on
string='legend(';
p_kstest_all_corrs = [];
edges=-0.2:0.004:1;
for i=1:NOF
    prop=histc(correlationsAll{i},edges);
    plot(edges,prop./sum(prop),'color',col{i})
    string=[string,'''Rec ',num2str(i),''','];
    if i>1
        [~,p] = kstest2(correlationsAll{i-1},correlationsAll{i});
        p_kstest_all_corrs = [p_kstest_all_corrs p];
    end
end
string=string(1:end-1);
string=[string,',''location'',''Northeast'')'];
eval(string)
xlabel('Correlation coefficient','fontsize',13)
ylabel('Histogram counts','fontsize',13)
title('HISTROGRAM','fontsize',14)
clear prop string edges i
[p_kruskal_corrs,anovatab_correlations] = kruskalwallis(correlations_pool,group_correlations,'off');
disp(['Nonparametric ANOVA for the Pearson correlation coefficients: ',num2str(p_kruskal_corrs)])
disp('     and Kolmogorov-Smirnov test between consecutive recordings:  ')
vpa(p_kstest_all_corrs,4)
disp(' ')
% now get the same data for channels that are active throughout the
% experiment:
clear col
[chankeyCommon,chankeyFilt]=GetSameChannels(chankeyAll);
rate_diff = [];
corr_diff = [];
% set colors:
cols = jet(NOF-1);
disp(' ')
disp(' Selecting common channels between all recordings ')
disp(' ')
disp(' ')
% now get the correlations and distances only for the channels of interest:
for i=1:NOF    
    % re-load the correlation data:
    eval(['load ',filenames{i}])
    distancesAll{i} = distances;
    correlationsAll{i} = correlations;
    % first find out which are kicked out:
    missing = [];
    prop = zeros(1,size(correlationsAll{i},1));
    for j=1:length(chankeyAll{i})
        prop = find((chankeyFilt{i}(:,1))==chankeyAll{i}(j,1));
        if isempty(prop)
            missing = [missing j];
        end
    end
    % now get rid of them in correlations and distances:
    correlationsAll{i}(missing,:)=[];
    correlationsAll{i}(:,missing)=[];
    distancesAll{i}(missing,:)=[];
    distancesAll{i}(:,missing)=[];
    % and pool together:
    distancesAll{i}=distancesAll{i}(:);
    correlationsAll{i}=correlationsAll{i}(:);
    correlationsAll{i}=correlationsAll{i}(distancesAll{i}~=0);
    distancesAll{i}=distancesAll{i}(distancesAll{i}~=0);    
    % get relative differences in rates and correlations in individual
    % channels (relative to the previous):
    if i>1
        rate_diff{i-1} = ( chankeyFilt{i}(:,4) - chankeyFilt{i-1}(:,4) ) ./ ( chankeyFilt{i}(:,4) + chankeyFilt{i-1}(:,4));
        corr_diff{i-1} =  ( correlationsAll{i} - correlationsAll{i-1} ) ./ ( abs(correlationsAll{i}) + abs(correlationsAll{i-1}) );
        col{i-1} = cols(i-1,:);
    end
end
% plots and stats of the changes in rates in individual channels:
figure(1)
subplot('position',[0.07 0.10 0.38 0.33])
PlotBoxplot(rate_diff,[],col);
xlabel('Recording','fontsize',13)
ylabel('Relative rate change','fontsize',13)
title('BOXPLOT','fontsize',14)
subplot('position',[0.56 0.10 0.4 0.35])
p_kstest_common_rates = [];
p_signtest_common_rates = [];
string=['legend('];
for i=1:NOF-1
    semilogx(chankeyFilt{i}(:,4)./rtime,rate_diff{i},'.','color',col{i})
    hold on
    [~,p] = kstest2(chankeyFilt{i+1}(:,4)./rtime,chankeyFilt{i}(:,4)./rtime);
    p_kstest_common_rates = [p_kstest_common_rates p];
    p = signtest(chankeyFilt{i+1}(:,4)./rtime,chankeyFilt{i}(:,4)./rtime);
    p_signtest_common_rates = [p_signtest_common_rates p];
    string=[string,'''Rec ',num2str(i+1),'vs',num2str(i),''','];
end
string=string(1:end-1);
string=[string,',''location'',''Northeast'')'];
eval(string)
xlabel('Log baseline firing rate','fontsize',13)
ylabel('Relative rate change','fontsize',13)
title('SCATTERPLOT','fontsize',14)
disp('Kolmogorov-Smirnov test for firing rates between consecutive recordings, after common channels were selected:  ')
vpa(p_kstest_common_rates,4)
disp('     and a paired signed rank test:  ')
vpa(p_signtest_common_rates,4)
disp(' ')
disp(' ')
% plots and stats of the changes in correlations in individual channels:
figure(2)
subplot('position',[0.07 0.10 0.38 0.33])
PlotBoxplot(corr_diff,[],col);
xlabel('Recording','fontsize',13)
ylabel('Relative correlation change','fontsize',13)
title('BOXPLOT','fontsize',14)
subplot('position',[0.56 0.10 0.4 0.35])
p_kstest_common_corrs = [];
p_signtest_common_corrs = [];
for i=1:NOF-1
    semilogx(abs(correlationsAll{i}),corr_diff{i},'.','color',col{i})
    hold on
    [~,p] = kstest2(correlationsAll{i+1},correlationsAll{i});
    p_kstest_common_corrs = [p_kstest_common_corrs p];
    p = signtest(correlationsAll{i+1},correlationsAll{i});
    p_signtest_common_corrs = [p_signtest_common_corrs p];
    
end
eval(string)
legend('location','Northwest')
xlabel('Log absolute baseline correlation','fontsize',13)
ylabel('Relative correlation change','fontsize',13)
title('SCATTERPLOT','fontsize',14)
disp(' ')
disp('Kolmogorov-Smirnov test for Pearson correlation coefficients between consecutive recordings, after common channels were selected:  ')
vpa(p_kstest_common_corrs,4)
disp('     and a paired signed rank test:  ')
vpa(p_signtest_common_corrs,4)
disp(' ')
%---------END----------%