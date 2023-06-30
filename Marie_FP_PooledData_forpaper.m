% Marie_FP_PooledData
% Marie Labouesse, marie.labouesse@gmail.com - Feb 2021

% 1- POOLED DATA
% load matlab spaces generated in: Marie_FP_IndivData_extraction
% puts individual trials into one big matrix (all mice together): PooledINDIV
% calculates the mean of individual trial trajectories aligned to specific event for each mouse: PooledAVE
% saves a workspace that can then be used for quantifications

% 2- POOLED GRAPHS
% generates average graphs curves (this used to be in the main script)
% generates heatmaps (1 average trace/mouse or all trials/all mice in one big heatmap - separated by power)

% 3- POOLED OVERALL CORRELATIONS: VALUES AND GRAPHS -
% calculates the mean of correlation values and graphs between FP data and behavior for each mouse across entire session
% generates average graphs and correlation value; commpare to shuffle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%% INITIALIZATIONS
close all; clear variables; clc;
set(0,'defaultfigurecolor',[1 1 1])

%% PARAMETERS 
show_plot = 1; % If 0, plots are not displayed
save_plot = 1; % If 0, plots are not saved
reanalysis = 0; % If 1, the code runs for sessions already analyzed. If 0, session excluded if analysed. This is based on the existence of the matlab space "IndividualData.mat"
overwrite = 0; % If reanalyzing data, 1 if you want to save them and overwrite results and 0 if not          
Collapse_sessions = 0; % If 0, analysis only of individual sessions. If 1, analysis of the group data (pooled)   

pooledtype = {'raw','baselinecorr'};

%% Define the path where the data is
PATH2DATA_0 = uigetdir('select folder'); %select the overall folder, where the Chrimson and mCherry are saved- we will import both genotypes
% virus = {'Chrimson','mCherry'}; 
virus = {'Chrimson'}; 
% virus = {'LED465'}; 
% virus = {'LED595'};

for v=1:length(virus)
    PATH2DATA.(virus{v}) = [PATH2DATA_0,'\',virus{v}];
    PATH2SAVEFOLDER.(virus{v}) = [PATH2DATA_0,'\',virus{v}];
    mice_list_virus.(virus{v}) = dir(PATH2DATA.(virus{v})); %all things in this folder
    mkdir([PATH2SAVEFOLDER.(virus{v}),'\pooled figures\']);
    mkdir([PATH2SAVEFOLDER.(virus{v}),'\pooled data\']);
end

                

%% IDENTIFY MICE TO ANALYZE 
for v=1:length(virus)
    for o = length(mice_list_virus.(virus{v})):-1:1
        if mice_list_virus.(virus{v})(o).isdir == 0  %remove non-folders
            mice_list_virus.(virus{v})(o) = [];
        else
            if  strcmp(mice_list_virus.(virus{v})(o).name,'data') == 1 || strcmp(mice_list_virus.(virus{v})(o).name,'figures') == 1 ...   %remove folders with data or figures or ., ..
                || contains(mice_list_virus.(virus{v})(o).name,'data') || contains(mice_list_virus.(virus{v})(o).name,'figures') ...
                || contains(mice_list_virus.(virus{v})(o).name,'results') || contains(mice_list_virus.(virus{v})(o).name,'turns') || contains(mice_list_virus.(virus{v})(o).name,'other') ...
                || strcmp(mice_list_virus.(virus{v})(o).name,'.') == 1 || strcmp(mice_list_virus.(virus{v})(o).name,'..') == 1
                mice_list_virus.(virus{v})(o) = [];
            end
        end
    end
    Nmice_virus{v} = length(mice_list_virus.(virus{v}));
%     AnimalID.(virus{v}) = [];
%     for nummice=1:length(mice_list_virus.(virus{v}))
%         AnimalID.(virus{v}){end+1} = ['ID',mice_list_virus.(virus{v})(nummice).name(end-3:end)];
%     end
end

              

%% Import individual workspaces and load the data into a pooled array
for v=1:length(virus)
    %AnimalIDs
    for nummice=1:length(mice_list_virus.(virus{v}))
        AnimalIDs(nummice) = {['ID',mice_list_virus.(virus{v})(nummice).name(end-3:end)]};
    end   
    
    for nummice=1:length(mice_list_virus.(virus{v}))
        % Define the path to the mouse and find the folders to analyze: 
        path2mouse = [PATH2DATA.(virus{v}),'\',mice_list_virus.(virus{v})(nummice).name,'\']; 
        % if there are several sessions inside the mouse folder
        sessions = dir(path2mouse);
        %remove non relevant folders
        for o = length(sessions):-1:1
            if sessions(o).isdir == 0  %remove non-folders
                sessions(o) = [];
            elseif strcmp(sessions(o).name,'data') == 1 || strcmp(sessions(o).name,'figures') == 1 || strcmp(sessions(o).name,'results') == 1 ...
                 || strcmp(sessions(o).name,'.') == 1 || strcmp(sessions(o).name,'..') == 1 || contains(sessions(o).name,'.csv') == 1 ...
                 || contains(sessions(o).name,'figures') == 1 || contains(sessions(o).name,'data') == 1 || contains(sessions(o).name,'other')  || contains(sessions(o).name,'BACKUP')
                 sessions(o) = [];
            end
        end
        if isempty(sessions)
            sessions(1).name = []; %to have an existent session
        end
        
        %SessionIDs
        for s=1:length(sessions)
            SessionIDs(s) = {['SessionNum',num2str(s)]};
            SessionNames(s) = {sessions(s).name};   
        end       
       
        % Create a folder to save the data for this mouse and define the path
        if exist([PATH2SAVEFOLDER.(virus{v}),'\',mice_list_virus.(virus{v})(nummice).name,'\'],'dir') == 0
            mkdir([PATH2SAVEFOLDER.(virus{v}),'\',mice_list_virus.(virus{v})(nummice).name,'\']);
        end
        path2save_mouse = [PATH2SAVEFOLDER.(virus{v}),'\',mice_list_virus.(virus{v})(nummice).name,'\'];
        
        %% Loop for all the sessions for the mouse: Sessions are experimental days you want to average together (replicates of each other)
        for s = 1:length(sessions)
            % Define the path to the session and create folder to save if needed:
            PATH2SESSION = [path2mouse,sessions(s).name];
            if exist([path2save_mouse,sessions(s).name],'dir') == 0
                mkdir([path2save_mouse,sessions(s).name])
            end
            if length(sessions) == 1
                PATH2SAVE = path2save_mouse;
            else
                PATH2SAVE = [path2save_mouse,sessions(s).name,'\'];
            end
            if exist([PATH2SAVE,'figures'],'dir') == 0
                mkdir([PATH2SAVE,'figures'])
            end
            if exist([PATH2SAVE,'results'],'dir') == 0
                mkdir([PATH2SAVE,'results'])
            end
        
            % Check if results are already saved for this session 
            done = exist([PATH2SAVE,'PooledAllMice.mat'],'file'); % .....
            if done == 0 || reanalysis == 1 %if not done or if wish to reanalyze
        
                % load the mouse/session matlab space
                load([PATH2SESSION,'\IndividualData.mat']);
       
                % Initialization of pooled structure with individual trials   %...

                if nummice == 1 && s == 1
                    for d=1:length(dFF_names)
                        for pow = 1:length(Opto_Powers)
                            for k = 1:size(datatype,2)
                                % ANIMAL ID ...................
                                PooledAnimalID.(virus{v}).(AnimalIDs{nummice}).(SessionIDs{s}) = {};
                                % INDIVIDUAL TRIALS OF EACH MICE FOR ALL MICE
                                PooledINDIV.(virus{v}).raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}) = ones(size(IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})))*nan;  
                                PooledINDIV.(virus{v}).baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}) = ones(length(mice_list_virus.(virus{v})),size(t_trials,2))*nan;  
                                % AVERAGE TRIAL OF EACH MICE (ALL SESSIONS) FOR ALL MICE, EACH MOUSE ON ONE LINE
                                PooledAVE.(virus{v}).raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}) = ones(length(mice_list_virus.(virus{v})),size(t_trials,2))*nan;  
                                PooledAVE.(virus{v}).baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}) = ones(length(mice_list_virus.(virus{v})),size(t_trials,2))*nan; 
                                % STREAMS                
                                fullstreams.(virus{v}).(AnimalIDs{nummice}).(SessionIDs{s}) = {};

                            end
                        end
                    end
                end
                
                % Add data one animal at a time
                % ANIMAL ID and SESSION ID
                PooledAnimalID.(virus{v}).(AnimalIDs{nummice}).(SessionIDs{s}) = SessionNames(s);
                PooledAnimalID_only.(virus{v}) = AnimalIDs;
                % STORE INDIVIDUAL TRIALS
                for d=1:length(dFF_names)
                    for pow = 1:length(Opto_Powers)
                        for k = 1:size(datatype,2)
                            PooledINDIV.(virus{v}).raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}) = IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}); 
                            PooledINDIV.(virus{v}).baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}) = IndivStim_data.baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}); 
                        end
                    end
                end                            
                % STORE AVERAGE
                for d=1:length(dFF_names)
                    for pow = 1:length(Opto_Powers)
                        for k = 1:size(datatype,2)
                            if length(sessions) == 1
                                PooledAVE.(virus{v}).raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})(nummice,:) = nanmean(IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}),1); 
                                PooledAVE.(virus{v}).baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})(nummice,:) = nanmean(IndivStim_data.baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}),1); 
                            elseif length(sessions) > 1
                            end
                        end
                    end
                end
                %STORE STREAMS
                fullstreams.(virus{v}).(AnimalIDs{nummice}).(SessionIDs{s}) = streams;
                

                
                
            end
        end
    end
end

%% Save pooled data for all viruses together
% STORE INDIVIDUAL and AVERAGE TRIALS
save([PATH2DATA_0,'\PooledAllMice.mat'],'PooledINDIV','PooledAVE','virus','PooledAnimalID','PooledAnimalID_only','dFF_names','datatype','Opto_Powers',...
    't_trials','stim_duration','TRANGE','BASELINE_WIN','Events2Align2','dt_ds','sampling_rate_ds','length_data','fullstreams','time_vect','pooledtype');




%% Create plot with nanmean and nanstd, Save graph

%%%%%%%%%%% GRAPH PARAMETERS (USED FOR AVERAGE GRAPHS)
% limits2plot.all = [-3 3]; %If empty, automatically adjusted to the data or whatever you specified.                
% limits2plot.(IdChannel{1}) = []; %If empty, automatically adjusted to the data or whatever you specified.       
% limits2plot.(IdChannel{2}) = []; %If empty, automatically adjusted to the data or whatever you specified. 
for d=1:length(dFF_names)
    for k=1:length(datatype)
        limits2plot.raw.(dFF_names{d}).(datatype{k}) = [-20 20]; %for all, then correct for speed
        limits2plot.baselinecorr.(dFF_names{d}).(datatype{k}) = [-20 20]; %for all, then correct for speed
    end
    if isfield(limits2plot.raw.(dFF_names{d}),'ZScoredFF')
        limits2plot.raw.(dFF_names{d}).ZScoredFF = [-4 4];
        limits2plot.baselinecorr.(dFF_names{d}).ZScoredFF = [-4 4];
    end
    if isfield(limits2plot.raw.(dFF_names{d}),'speed')
        limits2plot.raw.(dFF_names{d}).speed = [0 0.4];
        limits2plot.baselinecorr.(dFF_names{d}).speed = [-0.3 0.3];
    end
end

%plot
color2plot = {'b','g','m','r','c','k'};
TrialType_Number = length(Opto_Powers);

for v=1:length(virus)
    % for k = 1:size(Body_parts,2)
    for k = 1:size(datatype,2)
        for d=1:length(dFF_names)
            for p = 1:length(pooledtype)
                figure; clf; hold on
                if show_plot == 0
                    set(gcf,'visible','off')
                end
                for pow = 1:length(Opto_Powers)
                    tmp_avg_pool = nanmean(PooledAVE.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}),1); 
                    tmp_error_pool = nanstd(PooledAVE.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}),1,1)./...
                    sqrt(sum(~isnan(PooledAVE.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})(:,1))));
                    error_area(t_trials,tmp_avg_pool,tmp_error_pool,color2plot{pow},0.25); %t_trials is the time vector
                    max_avg(pow) = max(tmp_avg_pool);
                end
                % plot properties
                xline(0,'-k');
                xline(stim_duration,'-m');
                yline(0,'-.k');
                %     yline(max(max_avg),'-.r');
                xlabel('Time (s)');
                ylabel(datatype{k});
                if i == 4; % for speed
                    ylabel([datatype{k},' m/s']);
                end
                xlim([TRANGE(1) TRANGE(2)]);
%                 ylim(limits2plot.(pooledtype{p}).(dFF_names{d}).(datatype{k})); %see parameters setting at start                     

                PeakAve = num2str(mean(max_avg));
                PeakMax = num2str(max(max_avg));
                sgtitle([(pooledtype{p}),' ',datatype{k},', AvePeak: ',PeakAve,', MaxPeak: ',PeakMax],'Interpreter','none')

                %legend
                if TrialType_Number == 2
                % 2 powers
                legend([Opto_Powers{1}(2:end),'error'],[Opto_Powers{1}(2:end)],...
                        [Opto_Powers{2}(2:end),'error'],[Opto_Powers{2}(2:end)],...
                        'Opto On','Opto off','Location','best','NumColumns',1);      
                elseif TrialType_Number == 3
                % 3 powers
                legend([Opto_Powers{1}(2:end),'error'],[Opto_Powers{1}(2:end)],...
                    [Opto_Powers{2}(2:end),'error'],[Opto_Powers{2}(2:end)],...
                    [Opto_Powers{3}(2:end),'error'],[Opto_Powers{3}(2:end)],...
                    'Opto On','Opto off','Location','best','NumColumns',1);
                elseif TrialType_Number == 4
                % 4 powers
                legend([Opto_Powers{1}(2:end),'error'],[Opto_Powers{1}(2:end)],...
                    [Opto_Powers{2}(2:end),'error'],[Opto_Powers{2}(2:end)],...
                    [Opto_Powers{3}(2:end),'error'],[Opto_Powers{3}(2:end)],...
                    [Opto_Powers{4}(2:end),'error'],[Opto_Powers{4}(2:end)],...
                    'Opto On','Opto off','Location','best','NumColumns',1);
                        % 5 powers
                        legend([Opto_Powers{1}(2:end),'error'],[Opto_Powers{1}(2:end)],...
                            [Opto_Powers{2}(2:end),'error'],[Opto_Powers{2}(2:end)],...
                            [Opto_Powers{3}(2:end),'error'],[Opto_Powers{3}(2:end)],...
                            [Opto_Powers{4}(2:end),'error'],[Opto_Powers{4}(2:end)],...
                            [Opto_Powers{5}(2:end),'error'],[Opto_Powers{5}(2:end)],...
                            'Opto On','Opto off','Location','northwest','NumColumns',1)
                        elseif TrialType_Number == 6
                        % 6 powers
                        legend([Opto_Powers{1}(2:end),'error'],[Opto_Powers{1}(2:end)],...
                            [Opto_Powers{2}(2:end),'error'],[Opto_Powers{2}(2:end)],...
                            [Opto_Powers{3}(2:end),'error'],[Opto_Powers{3}(2:end)],...
                            [Opto_Powers{4}(2:end),'error'],[Opto_Powers{4}(2:end)],...
                            [Opto_Powers{5}(2:end),'error'],[Opto_Powers{5}(2:end)],...
                            [Opto_Powers{6}(2:end),'error'],[Opto_Powers{6}(2:end)],...
                            'Opto On','Opto off','Location','northwest','NumColumns',1)
                        end
                %saveplot or not
                if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                    saveas(gcf,[PATH2SAVEFOLDER.(virus{v}),'\pooled figures\AVE all powers optostim',pooledtype{p},' ',datatype{k},'.tif']);
                    saveas(gcf,[PATH2SAVEFOLDER.(virus{v}),'\pooled figures\AVE all powers optostim',pooledtype{p},' ',datatype{k},'.fig']);
                end
            end
        end
    end
    % end
end






%% Create a giant heatmaps - all mice included, using the AVE
orderedmap = 1; % 0 to plot in normal order, 1 if you want to plot in the order of deg inhibition by sorting on the minimal value between 3 and 10 sec (for 10 sec stim...)

wantoverallheatmap = 0;
if wantoverallheatmap == 1
    % OVERALL HEATMAP
    for v=1:length(virus)
        for p = 2; %1:length(pooledtype)
            for d=1; %:length(dFF_names)
                for k = 1; %:size(datatype,2)
                    if k == 1
                        Color_scale = [-2 2];
                    else
                        Color_scale = [-10 10];
                    end
                    figure;
                    if show_plot == 0
                        set(gcf,'visible','off')
                    end
                    for pow = 1:length(Opto_Powers)
                        subplot(length(Opto_Powers),1,pow)
                        % mean of individual trials for each mouse- into a heatmap that includes all mice
                        mean_array = PooledAVE.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k});
    %                     subplot(length(IdChannel)+1,1,o)
                        % if you want mice sorted by order in terms of deg inhibition
                        if orderedmap == 1
                            [~,order] = sort(abs(min(mean_array(:,t_trials > 2 & t_trials < stim_duration),[],2))); %if A is a matrix, then min(A,[],2) is a column vector containing the minimum value of each row
                            %B = sort(A) sorts the elements of A in ascending order + [B,I] = sort(___) also returns a collection of index vectors
                        else
                            order = [1:numrows]';
                        end
                        imagesc(t_trials,1,mean_array(order,:))
                        colormap('jet')
                        c = colorbar;
                        c.Label.String = 'Intensity (A.U.)';
                        if ~isempty(Color_scale)
                            c_limits = Color_scale;
                            caxis(c_limits)
                        end
                        hold on
                        xlim([t_trials(1) t_trials(end)])
                        xlabel('Time (s)')
                        ylabel('Animal')
                        set(gca,'YTick',1:size(mean_array,1));
                        title([Opto_Powers{pow}(2:end-2),' uW']);
                        box off
                        ylimits = get(gca,'YLim');
                        plot([0 0],ylimits,'k','LineWidth',2)
                        plot([stim_duration stim_duration],ylimits,'k','LineWidth',2,'LineStyle',':')

                        %     xline(0,'-k');
    %                     plot(round(t_Dipper(1,2)-t_Dipper(1,1)),ylimits,'m','LineWidth',2)
                    end
                    sgtitle([virus{v},', ',datatype{k},', ',pooledtype{p},', ',dFF_names{d}]);
                    if save_plot == 1
                        saveas(gcf,[PATH2SAVEFOLDER.(virus{v}),'\pooled figures\Heatmaps of dFF signals',pooledtype{p},' ',datatype{k},'.tif']);
                        saveas(gcf,[PATH2SAVEFOLDER.(virus{v}),'\pooled figures\Heatmaps of dFF signals',pooledtype{p},' ',datatype{k},'.fig']);
                    end
                end
            end
        end
    end
end

%% Individual trials- heatmaps
orderedmap =0; % 0 to plot in normal order, 1 if you want to plot in the order of deg inhibition by sorting on the minimal value between 3 and 10 sec (for 10 sec stim...)
for v=1:length(virus)
    for p = 2; %1:length(pooledtype)
        for d=1; %:length(dFF_names)
            for k = 1:size(datatype,2)
                if k == 1
                    Color_scale = [-1.1 1.4]; % Anymaze 0.2mW
%                     Color_scale = [-3 3]; % closed loop dFF
                    Color_scale = [-2 2]; % Anymaze 0.2mW

                else
                    Color_scale = [-10 10];
                end
                f= figure; f.Position = [100 100 600 550];
                if show_plot == 0
                    set(gcf,'visible','off')
                end
                for pow = 1:length(Opto_Powers)
                    subplot(length(Opto_Powers),1,pow);
                    alltrials_array = [];
                    for nummice=1:length(mice_list_virus.(virus{v}))
                        animals = PooledAnimalID_only.(virus{v});
                        for s = 1:length(sessions)
                            indivdata = PooledINDIV.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(animals{nummice}).(SessionIDs{s});
                            numrows = size(indivdata,1);
                            if isempty(alltrials_array)                               
                                alltrials_array(1:numrows,:) = indivdata;
                            else
                               numcurrentrows = size(alltrials_array,1);
                               alltrials_array(numcurrentrows+1:numcurrentrows+numrows,:) = indivdata;
                            end
                        end
                    end
                    % if you want each animal to have its trials sorted by order in terms of deg inhibition
                    if orderedmap == 1
                        [~,order] = sort(abs(min(alltrials_array(:,t_trials > 0.5 & t_trials < stim_duration),[],2))); %if A is a matrix, then min(A,[],2) is a column vector containing the minimum value of each row
                        %B = sort(A) sorts the elements of A in ascending order + [B,I] = sort(___) also returns a collection of index vectors
                    else
                        order = [1:size(alltrials_array,1)]';
                    end
                    %plot
                    imagesc(t_trials,1,alltrials_array(order,:))
                    colormap('redbluecmap')
                    c = colorbar;
                    c.Label.String = 'Zscore dFF';
                    set(gca,'FontSize',16,'FontName', 'Arial')  
                    ax = gca;
                    if ~isempty(Color_scale)
                        c_limits = Color_scale;
                        caxis(c_limits)
                    end
                    hold on
                    xlim([t_trials(1) t_trials(end)])
                    ylabel('Trial')
                    set(gca,'YTick',10:20:size(alltrials_array,1));
                    if pow==1; titleplan='0'; elseif pow==2; titleplan='0.2'; elseif pow==3; titleplan='2'; end;
                    title([titleplan,' mW'],'FontSize',14);
                    box off
                    ylimits = get(gca,'YLim');
                    plot([0 0],ylimits,'k','LineWidth',2)
                    plot([stim_duration stim_duration],ylimits,'k','LineWidth',2,'LineStyle',':')
                    clear alltrials_array
                    
                end
                xlabel('Time (s)')
                sgtitle([virus{v},', ',datatype{k},', ',pooledtype{p},', ',dFF_names{d}]);

                if save_plot == 1
                    saveas(f,[PATH2SAVEFOLDER.(virus{v}),'\pooled figures\Heatmaps of dFF signals',pooledtype{p},' ',datatype{k},'.tif']);
                    saveas(f,[PATH2SAVEFOLDER.(virus{v}),'\pooled figures\Heatmaps of dFF signals',pooledtype{p},' ',datatype{k},'.fig']);
                end
            end
        end
    end
 end
                    
                


%% Generate vectors every second

% draft to get cumulative angle every second
ix_select_during_if1515sec= [1527 1629 1731 1833 1934 2036 2138 2240 2341 2443 2545 2646 2748 2850 2952 3053]
ix_select_all = ix_select_during_if1515sec;
PooledAVE.Chrimson.baselinecorr.GPe.p0uW.angle(:,ix_select_all)
PooledAVE.Chrimson.baselinecorr.GPe.p200uW.angle(:,ix_select_all)


        