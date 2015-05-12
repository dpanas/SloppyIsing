%-------------------------------%
% script:   CompareFisher_Projections
%           Script for further analysis of Ising and FIM results, across
%           recordings from different days - looking at the evolution of
%           projections of parameters onto the eigenvectors of FIM;
%           IMPORTANT: projections are onto the eigenvectors of FIM from
%           the FIRST RECORDING
%
% dependancy: SetUpFigure, PlotWithShade, SetPanelAppearance, colormaps.mat
%
%
% DAP Oct 2014
% !!! no error control !!!
%-------------------------------%

% first clear the workspace:
clear

% variables are here:
%-------------------%
% and filenames (check not to overwrite!!!)
filenames{1} = './results_fisher/chip136_0_fisher_8n_new_prob10_bin5_filt.mat';
filenames{2} = './results_fisher/chip136_1_fisher_8n_new_prob10_bin5_filt.mat';
filenames{3} = './results_fisher/chip136_2_fisher_8n_new_prob10_bin5_filt.mat';
filenames{4} = './results_fisher/chip136_4_fisher_8n_new_prob10_bin5_filt.mat';
filenames{5} = './results_fisher/chip136_5_fisher_8n_new_prob10_bin5_filt.mat';
filenames{6} = './results_fisher/chip136_6_fisher_8n_new_prob10_bin5_filt.mat';
% optional:(default setting ,0, is comparison to first recording):
previous = 0;  
%-------------------%

% (below no more variables, just code)
%-------------------%
close all
addpath('./functions/')
load ./functions/colormaps.mat
NOF = length(filenames);
eval(['col = colorcell_',num2str(NOF),';'])
disp(' ')
disp(['Script comparing the projections of parameters onto eigenvectors of the Fisher Information Matrix.'])
disp(' ')
for i_file=1:NOF
    filename_fisher = filenames{i_file};
    eval(['load ',filename_fisher])
    disp(['Loaded Fisher results: ',filename_fisher])
    nos = length(parameter_matrix);
    non = length(parameter_matrix{1});
    % obtaining the eigenvectors in parameter space defined by FIM from
    % first recording for doing all projections:
    if i_file==1
        eigenvec_for_projections = eigenvectors;
    end
    % computing scalar projections on all eigenvectors:
    for i=1:nos
        % obtaining parameters as a vector (corresponding to eigenvectors):
        param = parameter_matrix{i}(1:non+1:end);
        for j=1:non
            param = [param parameter_matrix{i}(j,j+1:end)];  
        end
        for j=1:length(fisher_matrix{i})
            param_projection(i,j) = sum(param.*eigenvec_for_projections{i}(:,j)');
        end
        clear param
    end
    eval(['projections_',num2str(i_file),' = param_projection;'])
    
end
% and plot a figure for each pairwise comparison:
colorjet = jet(size(projections_1,2));
for i_file=2:NOF
    [fhandle,subplots]=SetUpFigure(2);
    eval(['subplot(1',num2str(NOF-1),num2str(i_file-1),')'])
    if previous==0
        eval(['proj_diff = abs( projections_',num2str(i_file),' - projections_1 );'])
    elseif previous==1
        eval(['proj_diff = abs( projections_',num2str(i_file),' - projections_',num2str(i_file-1),' );'])
    end
    eval(subplots{1})
    PlotWithShade(1:1:size(projections_1,2),mean(proj_diff),mean(proj_diff)-std(proj_diff),mean(proj_diff)+std(proj_diff),col{i_file},0.5);
    SetPanelAppearance('Average behavior','Eigenvector rank',[],[],'Projection change',[],[])
    eval(subplots{2})
    hold on
    for i=1:nos
        if previous==0
            eval(['scatter(projections_1(i,:),projections_',num2str(i_file),'(i,:),5,colorjet);']) 
        elseif previous==1
            eval(['scatter(projections_',num2str(i_file-1),'(i,:),projections_',num2str(i_file),'(i,:),5,colorjet);'])
        end
    end
    hold off
    eval(subplots{2})
    SetPanelAppearance('Scatterplot','baseline projection',[],[],'compared projection',[],[])
end
%---------END----------%