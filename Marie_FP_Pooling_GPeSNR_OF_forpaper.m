% Marie_FP_Pooling_GPeSNR_OF 
% Marie Labouesse, marie.labouesse@gmail.com - Jan 2022

% 1- POOLED DATA
% load matlab spaces generated in: Marie_FP_IndivData_extraction_GPeSNR_OF 
% need to call the folder in which all the animals are located (for now only works with 1 session/mouse)
% 1) puts entire traces of all animals (dFF, speed) data into one big matrix --> to run Pearson analyses with shuffle on all mice
% 2) trials
% puts the average of individual Measurements (maxima, etc) into one big matrix (all mice together): MeasurementsAVE --> for GraphPad
% calculates the mean of individual trial trajectories aligned to specific event for each mouse and puts into a matrix: PooledAVE --> for plotting

% 2- POOLED GRAPHS
% generates average graphs curves (this used to be in the main script)
% generates heatmaps (1 average trace/mouse or all trials/all mice in one big heatmap - separated by power)

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

color2plot = {[1.00,0.41,0.16],[0.72,0.27,1.00],[0.47,0.9,0.19],[1.0 0.1 0.3],[0 0 0.7],[0.5 0.5 0.5]}; % orange, purple, green,red, blue, grey 
color2plot_new = {[0.72,0.27,1.00],[0.52,0.27,1.00],[0.47,0.9,0.19],[0,0.7,0.19]}; % purple, purple-ish, green, green-ish



%% Define the path where the data is
PATH2DATA_0 = uigetdir('select folder'); %select the overall folder, in which there are folders where the 2 groups are saved
virus = {'GPeSNr'}; % one or two viruses/groups
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
                || contains(mice_list_virus.(virus{v})(o).name,'results') || contains(mice_list_virus.(virus{v})(o).name,'other') ...
                || strcmp(mice_list_virus.(virus{v})(o).name,'.') == 1 || strcmp(mice_list_virus.(virus{v})(o).name,'..') == 1
                mice_list_virus.(virus{v})(o) = [];
            end
        end
    end
    Nmice_virus{v} = length(mice_list_virus.(virus{v}));

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
                || strcmp(sessions(o).name,'anymaze') || strcmp(sessions(o).name,'representative figures') ...
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
        
        %% Loop for all the sessions for the mouse: Sessions are experimental days you want to average together 
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
        
               done = 0;

            if done == 0 || reanalysis == 1 %if not done or if wish to reanalyze
        
                %% Load the mouse/session matlab space
%                 load([PATH2SESSION,'\IndividualData_filtered.mat']);
                load([PATH2SESSION,'\IndividualData.mat']);
                
               
                
                
                
                %% Initialization of the 5 pooled structure with individual trials,assume only 1 session/mouse
                % 5 pooled structures: PooledStim_data, PooledStreams, PooledBodySpeed, PooledMeasurements/PooledMeasurements_con
                % All hold the average of the data (except the streams and speed where there is only 1 trial). 
                what2measure = {'maxima','minima','minimalatency','maximalatency'};
                what2measure2 = {'maxima','minima'};
                if nummice == 1 && s == 1
                    % ANIMAL ID 
                    PooledAnimalID = {};
                    for d=1:length(dFF_names)                       
                        % STREAMS
                        Pooledstreams.dFF.(dFF_names{d}) = ones(length(mice_list_virus.(virus{v})),length_data*2)*nan;
                        PooledBodySpeed = ones(length(mice_list_virus.(virus{v})),length_data*2)*nan;
                        % STIM DATA
                        for i=1:length(Events)
                            for r=1:length(raw_or_corr)
                                for q=1:length(what2measure)
                                    PooledStim_data.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r}) = ones(length(mice_list_virus.(virus{v})),size(t_trials,2))*nan;
                                    PooledStim_data.(Events{i}).body_speed.(raw_or_corr{r}) = ones(length(mice_list_virus.(virus{v})),size(t_trials,2))*nan;
                                    % MEASUREMENTS
                                    PooledMeasurements.(what2measure{q}).(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r}) = ...
                                        ones(length(mice_list_virus.(virus{v})),1)*nan; 
                                    PooledMeasurements.(what2measure{q}).(Events{i}).body_speed.(raw_or_corr{r}) = ...
                                        ones(length(mice_list_virus.(virus{v})),1)*nan;
                                end
                                for q=1:length(what2measure2)
                                    PooledMeasurements_con.(what2measure{q}).(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r}) = ...
                                        ones(length(mice_list_virus.(virus{v})),1)*nan;
                                    PooledMeasurements_con.(what2measure{q}).(Events{i}).body_speed.(raw_or_corr{r}) = ...
                                        ones(length(mice_list_virus.(virus{v})),1)*nan;    
                                end
                            end
                        end
                    end
                end
                
                %% Add data one animal at a time
                % ANIMAL ID
                PooledAnimalID{nummice} = AnimalIDs{nummice};
                for d=1:length(dFF_names)                       
                % STREAMS
                    Pooledstreams.dFF.(dFF_names{d})(nummice,1:length(streams.ZScoredFF.(dFF_names{d}))) = streams.ZScoredFF.(dFF_names{d}); % 
                    PooledBodySpeed(nummice,1:length(streams.ZScoredFF.(dFF_names{d}))) = BodySpeed;
                    % STIM DATA
                    for i=1:length(Events)
                        for r=1:length(raw_or_corr)
                            PooledStim_data.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(nummice,:) = nanmean(Stim_data.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r}),1);
                            PooledStim_data.(Events{i}).body_speed.(raw_or_corr{r})(nummice,:) = nanmean(Stim_data.(Events{i}).body_speed.(raw_or_corr{r}),1);
                            for q=1:length(what2measure)
                                % MEASUREMENTS
                                PooledMeasurements.(what2measure{q}).(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(nummice,1) = ...
                                    nanmean(Measurements.(what2measure{q}).(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r}),1);                                
                                PooledMeasurements.(what2measure{q}).(Events{i}).body_speed.(raw_or_corr{r})(nummice,1) = ...
                                    nanmean(Measurements.(what2measure{q}).(Events{i}).body_speed.(raw_or_corr{r}),1);
                            end
                            for q=1:length(what2measure2)
                                PooledMeasurements_con.(what2measure2{q}).(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(nummice,1) = ...
                                    nanmean(Measurements_con.(what2measure2{q}).(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r}),1);
                                PooledMeasurements_con.(what2measure2{q}).(Events{i}).body_speed.(raw_or_corr{r})(nummice,1) = ...
                                    nanmean(Measurements_con.(what2measure2{q}).(Events{i}).body_speed.(raw_or_corr{r}),1);
                            end
                        end
                    end
                end                
            end
        end
    end
end

            
                
%% Add a PooledMeasurements calculated on the average traces if you decide to. 
MeasurementsfromAverage = 1;
% allocate
ix11 = find(t_trials==0); % 0
periodend = 1.2; % 1.2 seconds
tp = find(t_trials-periodend>0);
ix21 = tp(1);
periodend = 1.5; % 1.5 seconds
tp = find(t_trials-periodend>0);
ix12 = tp(1);
periodend = 4.5; % 4.5 seconds
tp = find(t_trials-periodend>0);
ix22 = tp(1);
periodend = -0.8; % -0.8 seconds
tp = find(t_trials-periodend>0);
ix13 = tp(1);
periodend = 0.8; % 0.8 seconds
tp = find(t_trials-periodend>0);
ix23 = tp(1);
periodend = 2; % 2 seconds
tp = find(t_trials-periodend>0);
ix14 = tp(1);
periodend = 5; % 5 seconds
tp = find(t_trials-periodend>0);
ix24 = tp(1);
if MeasurementsfromAverage == 1
    for nummice=1:length(mice_list_virus.(virus{1}))
        for i=1:length(Events)
            for r=1:length(raw_or_corr)
                for d=1:length(dFF_names)                       
                    % MINIMA OR MAXIMA
                    if i==1 % onset
                        % maxima in dFF reached after mov onset in the later +1.5 to +4.5 period
                        [pks,locs] = findpeaks(PooledStim_data.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(nummice,ix12:ix22));
                        PooledMeasurements_AVE.maxima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(nummice) = ...
                            max(pks,[],2);
                        clear pks locs
                        % maxima in body speed reached after mov onset in the later +1.5 to +4.5 period                        
                        [pks,locs] = findpeaks(-PooledStim_data.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(nummice,ix11:ix21));
                        if isempty(pks) % if we can't find a minima then we take the value at 0 for the minima
                            pks=PooledStim_data.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(nummice,ix11);
                        end
                        PooledMeasurements_AVE.minima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(nummice) = ...
                            min(-pks,[],2);
                        clear pks locs
                        % maxima in body speed reached after mov onset in the later +1.5 to +4.5 period                        
                        [pks,locs] = findpeaks(PooledStim_data.(Events{i}).body_speed.(raw_or_corr{r})(nummice,ix12:ix22));
                        PooledMeasurements_AVE.maxima.(Events{i}).body_speed.(raw_or_corr{r})(nummice) = ...
                            max(pks,[],2);
                        clear pks locs
                        % minima in body speed reached after mov onset at time point 0                                        
                        PooledMeasurements_AVE.minima.(Events{i}).body_speed.(raw_or_corr{r})(nummice) = PooledStim_data.(Events{i}).body_speed.(raw_or_corr{r})(nummice,ix11); %value when x=0
                    else % offset
                        % minima in dFF reached after mov offset in the later +2 to +5 period                                                
                        [pks,locs] = findpeaks(-PooledStim_data.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(nummice,ix14:ix24));
                        if isempty(pks)  % then we take the minima
                            pks=min(PooledStim_data.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(nummice,ix14:ix24),[],2);
                        end
                        PooledMeasurements_AVE.minima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(nummice) = ...
                            min(-pks,[],2);
                        clear pks locs
                        % maxima in dFF reached after mov offset in the earlier -0.8 to +0.8 sec period                                                                                                
                        [pks,locs] = findpeaks(PooledStim_data.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(nummice,ix13:ix23));
                        PooledMeasurements_AVE.maxima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(nummice) = ...
                            max(pks,[],2);
                        clear pks locs
                        % minima in dFF reached after mov offset in the 1.5 to +4.5 period                                                                       
                        [pks,locs] = findpeaks(-PooledStim_data.(Events{i}).body_speed.(raw_or_corr{r})(nummice,ix12:ix22));
                        if isempty(pks)  % then we take the minima
                            pks=min(PooledStim_data.(Events{i}).body_speed.(raw_or_corr{r})(nummice,ix12:ix22),[],2);
                        end
                        PooledMeasurements_AVE.minima.(Events{i}).body_speed.(raw_or_corr{r})(nummice) = ...
                            min(-pks,[],2);
                        clear pks locs
                        % maxima in body speed reached after mov onset at time point 0                                                                
                        PooledMeasurements_AVE.maxima.(Events{i}).body_speed.(raw_or_corr{r})(nummice) = PooledStim_data.(Events{i}).body_speed.(raw_or_corr{r})(nummice,ix11); %value when x=0
                    end
                    % AMPLITUDE CHANGE
                    if i==1 % onset
                        PooledMeasurements_AVE.amplitudechange.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(nummice) = ...
                        PooledMeasurements_AVE.maxima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(nummice) - ...
                        PooledMeasurements_AVE.minima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(nummice);
                        PooledMeasurements_AVE.amplitudechange.(Events{i}).body_speed.(raw_or_corr{r})(nummice) = ...
                        PooledMeasurements_AVE.maxima.(Events{i}).body_speed.(raw_or_corr{r})(nummice) - ...
                        PooledMeasurements_AVE.minima.(Events{i}).body_speed.(raw_or_corr{r})(nummice);
                    else
                        PooledMeasurements_AVE.amplitudechange.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(nummice) = ...
                        PooledMeasurements_AVE.minima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(nummice) - ...
                        PooledMeasurements_AVE.maxima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(nummice);
                        PooledMeasurements_AVE.amplitudechange.(Events{i}).body_speed.(raw_or_corr{r})(nummice) = ...
                        PooledMeasurements_AVE.minima.(Events{i}).body_speed.(raw_or_corr{r})(nummice) - ...
                        PooledMeasurements_AVE.maxima.(Events{i}).body_speed.(raw_or_corr{r})(nummice);                        
                    end
                end
            end
        end
    end
end

 


%% Save pooled data for all viruses together
% STORE INDIVIDUAL and AVERAGE TRIALS
save([PATH2DATA_0,'\PooledAllMice.mat'],'virus','PooledAnimalID','dFF_names','PooledMeasurements',...
    'PooledMeasurements_con','Pooledstreams','PooledBodySpeed','PooledStim_data',...
    't_trials','TRANGE','BASELINE_WIN','dt_ds','sampling_rate_ds','length_data','time_vect');



%% SAVE FOLDER
% we only have 1 virus
PATH2SAVE = [PATH2SAVEFOLDER.(virus{v}),'\'];

% for v=1:length(virus)



%% TRIGGERED DATASETS NORMAL GRAPHS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% full plot dff, (first derivative) and speed
for r=1:length(raw_or_corr)
%             for r=2
%                 for d=1:length(dFF_names)
        figure; 
        for i=1:length(Events) 
            subplot(1,2,i);
            if show_plot == 0
                set(gcf,'visible','off')
            end
            hold on
            % dFF GPe
            d=1;
            yyaxis left
            tmp_avg = nanmean(PooledStim_data.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r}),1);
            tmp_error = nanstd(PooledStim_data.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r}),1,1)./...
                sqrt(sum(~isnan(PooledStim_data.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(:,1))));
            error_area(t_trials,tmp_avg,tmp_error,color2plot{2},0.25); % error_area(X, Y, barlength, color, alpha, varargin) to plot SEM
            % dFF SNr
            d=2;
%                         yyaxis left
            tmp_avg = nanmean(PooledStim_data.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r}),1);
            tmp_error = nanstd(PooledStim_data.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r}),1,1)./...
                sqrt(sum(~isnan(PooledStim_data.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(:,1))));
            error_area(t_trials,tmp_avg,tmp_error,color2plot{3},0.25); % error_area(X, Y, barlength, color, alpha, varargin) to plot SEM

                             
            if i==1 % mov onset
                if r==2 % baseline corr
%                     min_x=-5; max_x=5; min_y=-0.5; max_y=1.5; axis([min_x max_x min_y max_y]);
                    min_x=-5; max_x=5; min_y=-0.167; max_y=1; axis([min_x max_x min_y max_y]); 
                else % raw
                    min_x=-5; max_x=5; min_y=-0.5; max_y=0.7; axis([min_x max_x min_y max_y]);                                 
                end
            else % mov offset
                if r==2
%                     min_x=-5; max_x=5; min_y=-1.5; max_y=0.5; axis([min_x max_x min_y max_y]);
                    min_x=-5; max_x=5; min_y=-0.8; max_y=0.115; axis([min_x max_x min_y max_y]);                             
                else
                    min_x=-5; max_x=5; min_y=-0.2; max_y=0.8; axis([min_x max_x min_y max_y]);                                                             
                end
            end
            ylabel('Zscore dFF'); 

            % speed
            yyaxis right
            tmp_avg3 = nanmean(PooledStim_data.(Events{i}).body_speed.(raw_or_corr{r}),1);
            tmp_error3 = nanstd(PooledStim_data.(Events{i}).body_speed.(raw_or_corr{r}),1,1)./...
                sqrt(sum(~isnan(PooledStim_data.(Events{i}).body_speed.(raw_or_corr{r})(:,1))));
            error_area(t_trials,tmp_avg3,tmp_error3,color2plot{1},0.25);
            ylabel('Mouse speed (cm/s)'); 
            if i==1
                if r==2
%                     min_x=-5; max_x=5; min_y=-3; max_y=9; axis([min_x max_x min_y max_y]); 
                    min_x=-5; max_x=5; min_y=-1; max_y=6; axis([min_x max_x min_y max_y]); 
                else
                    min_x=-5; max_x=5; min_y=1; max_y=8; axis([min_x max_x min_y max_y]); 
                end
            else
               if r==2
%                     min_x=-5; max_x=5; min_y=-12; max_y=4; axis([min_x max_x min_y max_y]); 
                    min_x=-5; max_x=5; min_y=-7; max_y=1; axis([min_x max_x min_y max_y]); 
               else
                    min_x=-5; max_x=5; min_y=2; max_y=11; axis([min_x max_x min_y max_y]);                                
               end
            end
%                         max_avg(pow) = max(tmp_avg);
            xline(0,'-.k'); %other option: plot([0 0],limits2plot.dLight,'k')
%                         xline(stim_duration,'-m');
            yline(0,'-.k');
%                         yline(max(max_avg),'-.r');

%                         if i == 4; % for speed
%                             ylabel([datatype{i},' m/s']);
%                         end
            xlim([TRANGE(1) TRANGE(2)]);
%                         ylim(limits2plot.(raw_or_corr{r}).(dFF_names{d}).(datatype{i})); %see parameters setting at start                     
%                         PeakAve = num2str(mean(max_avg));
%                         PeakMax = num2str(max(max_avg));
            sgtitle([raw_or_corr{r}],'Interpreter','none')
%                         sgtitle([raw_or_corr{r},' ',datatype{i},' ',dFF_names{d},', AvePeak: ',PeakAve,', MaxPeak: ',PeakMax],'Interpreter','none')

            % legend
            if i==1
            legend('dFF GPe',' ','dFF SNr',' ','Mouse speed',' ','Location','Northwest','NumColumns',1)      
            legend boxoff
            end
%                         legend('dFF error',' ','dFF first derivative', ' ','Mouse speed',' ','Location','best','NumColumns',1)      

            % title
            title(['movement ',Events{i}(5:end)],'Fontsize',12);

            % size 
            set(gca,'FontSize',16)  
            ax = gca;
            ax.YAxis(1).Color = 'k';
            ax.YAxis(2).Color = 'k';

            %saveplot or not
            if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                saveas(gcf,[PATH2SAVE,'figures\AVE movement dff ',raw_or_corr{r},' ',datatype{i},' ',dFF_names{d},'.tif'])
                saveas(gcf,[PATH2SAVE,'figures\AVE movement dff ',raw_or_corr{r},' ',datatype{i},' ',dFF_names{d},'.fig'])
            end
        end
%                 end
end %for d=1:length(dFFnames)



%% POOLED CORRELATION ANALYSES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Trim so that all Pooledstreams and PooledBodySpeed are the same duration, as well as time_vect   

idx_end_all = NaN(size(Pooledstreams.dFF.(dFF_names{d}),1),1);
for i=1:size(Pooledstreams.dFF.(dFF_names{d}),1)
   [row, col] = find(isnan(Pooledstreams.dFF.(dFF_names{d})(i,:))); 
   idx_end_all(i,1) = col(1);
end
idx_end = min(idx_end_all)

for d=1:length(dFF_names)                
    Pooledstreams.dFF.(dFF_names{d})(:,idx_end:end) = [];
    PooledBodySpeed(:,idx_end:end) = [];
    time_vect(:,idx_end:end) = [];
end




 %% Lag Analysis dFF GPe and SNr 
% set parameters
PATH2SAVE1 = PATH2SAVE;
Max_lag = 10; % Max time lag in seconds to test for correlation
corr_type = 'Pearson'; % Pearson or Spearman
save_name = 'Pearson';
save_plot = 0;
show_plot = 1; %set(gcf,'visible','off')
loop1 = 0;
calc_plot = 1;
show_laginfo = 1;
nummice = 1;
sample1 = Pooledstreams.dFF.(dFF_names{1})(nummice,:); % 
sample2 = Pooledstreams.dFF.(dFF_names{2})(nummice,:);  

% first lag2plot to get the lag2plot vector
[lag2plot,Ovrl_corr] = ovrl_corr_calculation_loop(sample1,sample2,time_vect,...
    Max_lag,corr_type,calc_plot,show_plot,save_plot,PATH2SAVE1,save_name,loop1,show_laginfo);

% pooled correlation: run once per animal, keep the data, then generate average+SEM plot
Pearson_output_GPESNR_real_max = NaN(length(mice_list_virus.(virus{1})),1); % column 1
Pearson_output_GPESNR_real_lag = NaN(length(mice_list_virus.(virus{1})),1); % column 1
lengthdataset = NaN(length(mice_list_virus.(virus{1})),1); % column 1
Pearson_output_GPESNR_real = NaN(length(mice_list_virus.(virus{1})),length(lag2plot)); % column 1
loop1 = 1;
show_laginfo = 0;
figure; clf; 
for nummice=1:length(mice_list_virus.(virus{1})) % running the pearson r correlation one animal at a time
    sample1 = Pooledstreams.dFF.(dFF_names{1})(nummice,:); % 
    sample2 = Pooledstreams.dFF.(dFF_names{2})(nummice,:); % 
    [lag2plot,Ovrl_corr] = ovrl_corr_calculation_loop(sample1,sample2,time_vect,...
    Max_lag,corr_type,calc_plot,show_plot,save_plot,PATH2SAVE1,save_name,loop1,show_laginfo);
%     plot(lag2plot,Ovrl_corr.r,'k'); hold on;
    Pearson_output_GPESNR_real(nummice,:) = Ovrl_corr.r;
    Pearson_output_GPESNR_real_max(nummice) = max(Ovrl_corr.r,[],2);
    Pearson_output_GPESNR_real_lag(nummice) = Ovrl_corr.positive.lag.val;
    lengthdataset(nummice)=length(sample1);
end

 
% VERSION 2: FULL SHUFFLE, incl PHASE SHUFFLE
% Then generate the shuffled population and use the last graph: 
numbershuffle_peranimal = 5; % if 125 times per mice, 8 mice, then its 1000 shuffles; CAREFUL TO TEST OUT THE CODE AT LOW NUMBER EG 5 TO NOT CRASH
minlengthdataset = min(lengthdataset,[],1);
Pearson_output_GPESNR_shuffle_max = NaN(length(mice_list_virus.(virus{1})),1); % column 1
Pearson_output_GPESNR_shuffle = NaN(length(mice_list_virus.(virus{1}))*numbershuffle_peranimal,length(lag2plot)); % column 1
save_plot = 0;
show_plot = 0; %set(gcf,'visible','off')
loop1 = 0;
calc_plot = 0;
show_laginfo = 0;
numbershuffletotal = numbershuffle_peranimal * length(mice_list_virus.(virus{1}));
for nummice=1:length(mice_list_virus.(virus{1})) % running the pearson r correlation one animal at a time
    w = waitbar(0, 'Preparing test...', 'Name', ['Shuffle, mouse ',num2str(nummice)]); 
    waitbar(0, w, sprintf('Shuffle 1 of %d', numbershuffle_peranimal), 'Name', ['Shuffle, mouse ',num2str(nummice)]); 
    for u=1:numbershuffle_peranimal
        waitbar(u/numbershuffle_peranimal, w, sprintf('Permutation %d of %d', u, numbershuffle_peranimal)); 
        
          
        % option 3: phase shuffle
        [sample1] = PhaseShuffle(Pooledstreams.dFF.(dFF_names{1})(nummice,:));
        index2 = [1:length(sample1)];
        
        % generate new samples
        sample2 = Pooledstreams.dFF.(dFF_names{2})(nummice,index2); %
        clear shuffleindex 
        [lag2plot,Ovrl_corr] = ovrl_corr_calculation_loop(sample1,sample2,time_vect,...
        Max_lag,corr_type,calc_plot,show_plot,save_plot,PATH2SAVE1,save_name,loop1,show_laginfo);
    %     plot(lag2plot,Ovrl_corr.r,'k'); hold on;
        if nummice==1
            rowindex2go = u;
        else
            rowindex2go = ((nummice-1) * numbershuffle_peranimal) + u;
        end
        Pearson_output_GPESNR_shuffle(rowindex2go,:) = Ovrl_corr.r;
        Pearson_output_GPESNR_shuffle_max(rowindex2go,1) = max(Ovrl_corr.r,[],2);
    end
    delete(w); 
end



% Getting probability of finding observed difference from random permutations (but the real stats are done later in GraphdPAD for the paper)
pvalue_realshuffle_GPeSNr = (length(find(mean(Pearson_output_GPESNR_real_max) < Pearson_output_GPESNR_shuffle_max))+1) / ((numbershuffle_peranimal*nummice)+1)


 
% Full plot
% first the error area for REAL data
figure; tmp_avg_GPeSNr = nanmean(Pearson_output_GPESNR_real,1);   
tmp_error_GPeSNr = nanstd(Pearson_output_GPESNR_real,1,1)./...
    sqrt(sum(~isnan(Pearson_output_GPESNR_real)));
error_area([lag2plot],tmp_avg_GPeSNr,tmp_error_GPeSNr,'r',0.25); % error_area(X, Y, barlength, color, alpha, varargin) to plot SEM  

% then the error area for SHUFFLE data; or choose to plot the 95% interval below
tmp_avg_GPeSNr_shuffle = nanmean(Pearson_output_GPESNR_shuffle,1);
tmp_error_GPeSNr_shuffle = nanstd(Pearson_output_GPESNR_shuffle,1,1)./...
    sqrt(sum(~isnan(Pearson_output_GPESNR_shuffle)));
% gcf; error_area([lag2plot],tmp_avg_GPeSNr_shuffle,tmp_error_GPeSNr_shuffle,'b',0.25); % error_area(X, Y, barlength, color, alpha, varargin) to plot SEM  


% or plot the 95% percentile, if normally distributed; Student's t distribution
CI95 = NaN(2,size(Pearson_output_GPESNR_shuffle,2));
for g=1:size(Pearson_output_GPESNR_shuffle,2)
    x = Pearson_output_GPESNR_shuffle(:,g)';
    ts = tinv([0.025 0.975],length(x)-1);      % T score
    SEM =  nanstd(x,1,2)./sqrt(sum(~isnan(x)));
    mean_x = nanmean(x,2);
    CI95(:,g) = mean(x) + ts*SEM;                      % Confidence Intervals
end
plot(lag2plot,nanmean(Pearson_output_GPESNR_shuffle,1),'b')
patch([lag2plot, fliplr(lag2plot)], [CI95(1,:) fliplr(CI95(2,:))], 'b' , 'EdgeColor' , 'none' , 'FaceAlpha' ,0.25)

% plot(lag2plot,Ovrl_corr.r,'k')
axis([[lag2plot(1) lag2plot(end)] -0.2 1]); 
hold on
plot([0 0],[-0.2 1],'--k')

xlabel('Lag (s)')
ylabel('Pearson r coefficient')
title('Correlation GPe/SNr')
box off


txt = ['\leftarrow + lag = ',num2str(nanmean(Pearson_output_GPESNR_real_lag),1),' s'];
text(nanmean(Pearson_output_GPESNR_real_lag),nanmean(Pearson_output_GPESNR_real_max)+0.07,txt,'FontSize',16)


txt2 = ['real vs shuffle: p = ',num2str(pvalue_realshuffle_GPeSNr,2)];
text(-5.5,-0.12,txt2,'Color','b','FontSize',16)

% Pearson level
r_ave = mean(Pearson_output_GPESNR_real_max);
yline(r_ave,'-.r','LineWidth',1)

% size 
set(gca,'FontSize',20,'FontName', 'Arial')  
ax = gca;
 
if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
    saveas(gcf,[PATH2SAVE,'figures\dFF Pearson corr GPe SNr.fig'])
    saveas(gcf,[PATH2SAVE,'figures\dFF Pearson corr GPe SNr.tif'])
end
% allcorr(z)=max(Ovrl_corr.r,[],2);



%% Correlations between dFF and Behavioral data, with smooth windows                      
        smoothwindows = [0.1 0.2 0.5 1 2 5 10 20 50]; % in seconds
for d=1:length(dFF_names)    
    allcorr.(dFF_names{d}) = NaN(length(mice_list_virus.(virus{1})),length(smoothwindows));
end
for nummice=1:length(mice_list_virus.(virus{1})) % running the pearson r correlation one animal at a time
    for d=1:length(dFF_names)                
        % Lag Analysis dFF and speed at different speed smoothing 
        smoothwindows = [0.1 0.2 0.5 1 2 5 10 20 50]; % in seconds
        for z= 1:length(smoothwindows)
            smooth_window = smoothwindows(z)./dt_ds; % first element here is in seconds
            smooth_window_sec = sprintf('%.0f',smooth_window*dt_ds)
            PooledBodySpeed_smooth = smooth(PooledBodySpeed(nummice,:),smooth_window);

            Max_lag = 10; % Max time lag in seconds to test for correlation
            corr_type = 'Pearson'; % Pearson or Spearman
            save_name = 'Pearson';
            sample1 = Pooledstreams.dFF.(dFF_names{d})(nummice,:); % 
            sample2 = PooledBodySpeed_smooth'; % take transpose if needed
            if show_plot == 0
               set(gcf,'visible','off')
            end
            calc_plot=1;
            loop0=0;
            show_laginfo = 1;
            [lag2plot,Ovrl_corr] = ovrl_corr_calculation_loop...
                (sample1,sample2,time_vect,Max_lag,corr_type,calc_plot,show_plot,save_plot,PATH2SAVE,save_name,loop0,show_laginfo);
            if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                saveas(gcf,[PATH2SAVE,'figures\dFF speed Pearson corr ',smooth_window_sec,'s smooth ',dFF_names{d},'.tif'])
                saveas(gcf,[PATH2SAVE,'figures\dFF speed Pearson corr ',smooth_window_sec,'s smooth',dFF_names{d},'.tif'])
            end
            allcorr.(dFF_names{d})(nummice,z)=max(Ovrl_corr.r,[],2);
        end      
    end
end

% average plot
color2plot_final = {[1.0 0.1 0.3],[0.72,0.27,1.00]}; % orange, purple, green,red, blue, grey 
f=figure; hold on
for d=1:length(dFF_names)    
    tmp_avg = nanmean(allcorr.(dFF_names{d}),1);
    tmp_error = nanstd(allcorr.(dFF_names{d}),1,1)./...
        sqrt(sum(~isnan(allcorr.(dFF_names{d}))));
    error_area([smoothwindows],tmp_avg,tmp_error,color2plot_final{d},0.25); % error_area(X, Y, barlength, color, alpha, varargin) to plot SEM  
    
    txt = [dFF_names{d},': max corr = ',num2str(max(tmp_avg,[],2),2),', at ',num2str(smoothwindows(find(tmp_avg==max(tmp_avg,[],2))),2),' sec bins'];
    %                     text(80,80,txt,'HorizontalAlignment','right')
    if  d==1; h=0.24; else; h=0.29; end
    text(1,h,txt,'FontSize',16,'Color',color2plot_final{d})
    scatter(smoothwindows(find(tmp_avg==max(tmp_avg,[],2))),max(tmp_avg,[],2),'LineWidth',1,'MarkerEdgeColor','k')
end
min_x=-1; max_x=50; min_y=0.2; max_y=0.7; axis([min_x max_x min_y max_y]); 
xlabel('Speed smoothing bin size (sec)'); ylabel('Pearson r'); 
title('dFF and speed correlation');
% size 
set(gca,'FontSize',20,'FontName', 'Arial')  
ax = gca;

if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
        saveas(gcf,[PATH2SAVE,'figures\dFF speed Pearson corr overall.fig'])
        saveas(gcf,[PATH2SAVE,'figures\dFF speed Pearson corr overall.tif'])
end


%% Correlations dFF speed, with shuffle

PATH2SAVE1 = PATH2SAVE;
Max_lag = 10; % Max time lag in seconds to test for correlation
corr_type = 'Pearson'; % Pearson or Spearman
save_name = 'Pearson';
save_plot = 0;
show_plot = 1; %set(gcf,'visible','off')
loop1 = 0;
calc_plot = 1;
show_laginfo = 1;
nummice = 1;
sample1 = Pooledstreams.dFF.(dFF_names{1})(nummice,:); % 
smoothwindows = 2; % in seconds
smooth_window = smoothwindows./dt_ds; % first element here is in seconds
sample2 = smooth(PooledBodySpeed(nummice,:),smooth_window)';

% first lag2plot to get the lag2plot vector
[lag2plot,Ovrl_corr] = ovrl_corr_calculation_loop(sample1,sample2,time_vect,...
    Max_lag,corr_type,calc_plot,show_plot,save_plot,PATH2SAVE1,save_name,loop1,show_laginfo);

% pooled correlation: run once per animal, keep the data, then generate average+SEM plot
clear Pearson_output_real_max Pearson_output_real_lag Pearson_output_real
Pearson_output_real_max.(dFF_names{d}) = NaN(length(mice_list_virus.(virus{1})),1); % column 1
Pearson_output_real_lag.(dFF_names{d}) = NaN(length(mice_list_virus.(virus{1})),1); % column 1
lengthdataset = NaN(length(mice_list_virus.(virus{1})),1); % column 1
Pearson_output_real.(dFF_names{d}) = NaN(length(mice_list_virus.(virus{1})),length(lag2plot)); % column 1
loop1 = 1;
show_laginfo = 0;
for d=1:length(dFF_names)    
figure; clf; 
    for nummice=1:length(mice_list_virus.(virus{1})) % running the pearson r correlation one animal at a time
        sample1 = Pooledstreams.dFF.(dFF_names{d})(nummice,:); % 
        sample2 = smooth(PooledBodySpeed(nummice,:),smooth_window)';
        [lag2plot,Ovrl_corr] = ovrl_corr_calculation_loop(sample1,sample2,time_vect,...
        Max_lag,corr_type,calc_plot,show_plot,save_plot,PATH2SAVE1,save_name,loop1,show_laginfo);
        Pearson_output_real.(dFF_names{d})(nummice,:) = Ovrl_corr.r;
        Pearson_output_real_max.(dFF_names{d})(nummice) = max(Ovrl_corr.r,[],2);
  
        lengthdataset(nummice)=length(sample1);
    end
end


% VERSION 2: PHASE SHUFFLE
% Then generate the shuffled population and use the last graph: 
clear Pearson_output_shuffle_max Pearson_output_shuffle_lag Pearson_output_shuffle
numbershuffle_peranimal = 2; % if 125 times per mice, 8 mice, then its 1000 shuffles CAREFUL TO TEST OUT THE CODE, RUN THIS AT LOW NUMBER EG 5 TO NOT CRASH
minlengthdataset = min(lengthdataset,[],1);
Pearson_output_shuffle_max.(dFF_names{d}) = NaN(length(mice_list_virus.(virus{1})),1); % column 1
Pearson_output_shuffle.(dFF_names{d}) = NaN(length(mice_list_virus.(virus{1}))*numbershuffle_peranimal,length(lag2plot)); % column 1
save_plot = 0;
show_plot = 0; %set(gcf,'visible','off')
loop1 = 0;
calc_plot = 0;
show_laginfo = 0; 
numbershuffletotal = numbershuffle_peranimal * length(mice_list_virus.(virus{1}));
for d=1:length(dFF_names)    
    for nummice=1:length(mice_list_virus.(virus{1})) % running the pearson r correlation one animal at a time
        w = waitbar(0, 'Preparing test...', 'Name', ['Shuffle, mouse ',num2str(nummice)]); 
        waitbar(0, w, sprintf('Shuffle 1 of %d', numbershuffle_peranimal), 'Name', ['Shuffle, mouse ',num2str(nummice)]); 
        for u=1:numbershuffle_peranimal
            waitbar(u/numbershuffle_peranimal, w, sprintf('Permutation %d of %d', u, numbershuffle_peranimal)); 

            % option 3: phase shuffle
            [sample1] = PhaseShuffle(Pooledstreams.dFF.(dFF_names{d})(nummice,:));
            index2 = [1:length(sample1)];

            % generate new samples
            sample2 = smooth(PooledBodySpeed(nummice,index2),smooth_window)';

            clear shuffleindex 
            [lag2plot,Ovrl_corr] = ovrl_corr_calculation_loop(sample1,sample2,time_vect,...
            Max_lag,corr_type,calc_plot,show_plot,save_plot,PATH2SAVE1,save_name,loop1,show_laginfo);
        %     plot(lag2plot,Ovrl_corr.r,'k'); hold on;
            if nummice==1
                rowindex2go = u;
            else
                rowindex2go = ((nummice-1) * numbershuffle_peranimal) + u;
            end
            Pearson_output_shuffle.(dFF_names{d})(rowindex2go,:) = Ovrl_corr.r;
            Pearson_output_shuffle_max.(dFF_names{d})(rowindex2go,1) = max(Ovrl_corr.r,[],2);
        end
        delete(w); 
    end
end


% Getting probability of finding observed difference from random permutations
clear pvalue_realshuffle
for d=1:length(dFF_names)    
    pvalue_realshuffle.(dFF_names{d}) = (length(find(mean(Pearson_output_real_max.(dFF_names{d})) < Pearson_output_shuffle_max.(dFF_names{d})))+1) / ((numbershuffle_peranimal*nummice)+1)
end

 
% Full plot
% first the error area for REAL data
figure; hold on;
for d=1:length(dFF_names)    
    tmp_avg_dffspeed.(dFF_names{d}) = nanmean(Pearson_output_real.(dFF_names{d}),1);   
    tmp_error_dffspeed.(dFF_names{d}) = nanstd(Pearson_output_real.(dFF_names{d}),1,1)./...
        sqrt(sum(~isnan(Pearson_output_real.(dFF_names{d}))));
    error_area([lag2plot],tmp_avg_dffspeed.(dFF_names{d}),tmp_error_dffspeed.(dFF_names{d}),color2plot_final{d},0.25); % error_area(X, Y, barlength, color, alpha, varargin) to plot SEM  

% then the error area for SHUFFLE data; or choose to plot the 95% interval below
    tmp_avg_shuffle.(dFF_names{d}) = nanmean(Pearson_output_shuffle.(dFF_names{d}),1);
    tmp_error_shuffle.(dFF_names{d}) = nanstd(Pearson_output_shuffle.(dFF_names{d}),1,1)./...
        sqrt(sum(~isnan(Pearson_output_shuffle.(dFF_names{d}))));
  
    
    % or plot the 95% percentile, if normally distributed; Student's t distribution
    CI95 = NaN(2,size(Pearson_output_shuffle.(dFF_names{d}),2));
    for g=1:size(Pearson_output_shuffle.(dFF_names{d}),2)
        x = Pearson_output_shuffle.(dFF_names{d})(:,g)';
        ts = tinv([0.025 0.975],length(x)-1);      % T score
        SEM =  nanstd(x,1,2)./sqrt(sum(~isnan(x)));
        mean_x = nanmean(x,2);
        CI95(:,g) = mean(x) + ts*SEM;                      % Confidence Intervals
    end
    color2plot_shuffle = {[0 0 0.6],[0 0 0.95]}; % 

    plot(lag2plot,nanmean(Pearson_output_shuffle.(dFF_names{d}),1),'b')
    patch([lag2plot, fliplr(lag2plot)], [CI95(1,:) fliplr(CI95(2,:))], color2plot_shuffle{d} , 'EdgeColor' , 'none' , 'FaceAlpha' ,0.25)

    % plot(lag2plot,Ovrl_corr.r,'k')
    axis([[lag2plot(1) lag2plot(end)] -0.2 0.8]); 
    hold on
    plot([0 0],[-0.2 1],'--k')
  
    xlabel('Lag (s)')
    ylabel('Pearson r coefficient')
    title('Correlation dFF/mouse speed')
    box off
 
    % Pearson level
    r_ave = mean(Pearson_output_real_max.(dFF_names{d}));
    yline(r_ave,'-.','Color',color2plot_final{d},'LineWidth',1)

    % size 
    set(gca,'FontSize',20,'FontName', 'Arial')  
    ax = gca;
    
txt = ['\leftarrow + lag = ',num2str(nanmean(Pearson_output_real_lag.(dFF_names{d})),2),' s'];
if d == 1; o = 0.05; else o = -0.05; end;
text(nanmean(Pearson_output_real_lag.(dFF_names{d}))+0.7,nanmean(Pearson_output_real_max.(dFF_names{d}))+o,txt,'FontSize',16,'Color',color2plot_final{d})

end

txt2 = ['real vs shuffle: p = ',num2str(pvalue_realshuffle.(dFF_names{1}),2),' / ',num2str(pvalue_realshuffle.(dFF_names{2}),2)];
text(-7.5,-0.12,txt2,'Color','b','FontSize',16) 
    
save_plot = 1;
if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
    saveas(gcf,[PATH2SAVE,'figures\dFF Pearson corr GPe SNr and speed.fig'])
    saveas(gcf,[PATH2SAVE,'figures\dFF Pearson corr GPe SNr and speed.tif'])
end
% allcorr(z)=max(Ovrl_corr.r,[],2);



%% Save pooled data for all viruses together
% STORE INDIVIDUAL and AVERAGE TRIALS
save([PATH2DATA_0,'\PooledAllMice_forgraphpad_notfiltered_GPeSNr_butnotdFFspeed.mat'],'virus','PooledAnimalID','dFF_names','PooledMeasurements',...
    'PooledMeasurements_con','Pooledstreams','PooledBodySpeed','PooledStim_data',...
    't_trials','TRANGE','BASELINE_WIN','dt_ds','sampling_rate_ds','length_data','time_vect',...
    'allcorr','Pearson_output_real','Pearson_output_real_lag','pvalue_realshuffle',...
    'Pearson_output_real_max','Pearson_output_shuffle','tmp_avg_shuffle','tmp_error_shuffle',...
    'Pearson_output_GPESNR_real','Pearson_output_GPESNR_real_max','Pearson_output_GPESNR_real_lag',...
    'Pearson_output_GPESNR_shuffle','Pearson_output_GPESNR_shuffle_max','tmp_avg_GPeSNr','tmp_error_GPeSNr',...
    'tmp_avg_GPeSNr_shuffle','tmp_error_GPeSNr_shuffle','tmp_avg_dffspeed','tmp_error_dffspeed',...
    'pvalue_realshuffle_GPeSNr','Ovrl_corr_dffspeed');




        