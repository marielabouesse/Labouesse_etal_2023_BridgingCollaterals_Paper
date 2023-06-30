% SyGCAMP6s_pool

% Marie Labouesse, marie.labouesse@gmail.com - June 2021

% 1- POOLED DATA
% load matlab spaces generated in: Marie_FP_IndivData_extraction
% puts individual trials into one big matrix (all mice together): PooledINDIV
% calculates the mean of individual trial trajectories aligned to specific event for each mouse: PooledAVE
% saves a workspace that can then be used for quantifications

% 2- POOLED GRAPHS
% generates average graphs curves (this used to be in the main script)
% generates heatmaps (1 average trace/mouse or all trials/all mice in one big heatmap)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%% INITIALIZATIONS
close all; clear variables; clc;
set(0,'defaultfigurecolor',[1 1 1])

%% PARAMETERS 
show_plot = 1; % If 0, plots are not displayed
save_plot = 1; % If 0, plots are not saved
reanalysis = 0; % If 1, the code runs for sessions already analyzed. If 0, session excluded if analysed. This is based on the existence of the matlab space "IndividualData.mat"
overwrite = 0; % If reanalyzing data, 1 if you want to save them and overwrite results and 0 if not          
Collapse_sessions = 0; % If 0, analysis only of individual sessions. 

%% Number of sessions
numberofsessions = 5;

%% Define the path where the data is
PATH2DATA = uigetdir('select folder');      
PATH2SAVEFOLDER = PATH2DATA;
mice_list_virus = dir(PATH2DATA); %all things in this folder
mkdir(PATH2SAVEFOLDER,'\pooled figures\');
mkdir(PATH2SAVEFOLDER,'\pooled data\');

                

%% IDENTIFY MICE TO ANALYZE 
for o = length(mice_list_virus):-1:1
    if mice_list_virus(o).isdir == 0  %remove non-folders
        mice_list_virus(o) = [];
    else
        if  strcmp(mice_list_virus(o).name,'data') == 1 || strcmp(mice_list_virus(o).name,'figures') == 1 ...   %remove folders with data or figures or ., ..
            || contains(mice_list_virus(o).name,'data') || contains(mice_list_virus(o).name,'figures') ...
            || contains(mice_list_virus(o).name,'results') || contains(mice_list_virus(o).name,'FP-Rotarod') ...
            || strcmp(mice_list_virus(o).name,'.') == 1 || strcmp(mice_list_virus(o).name,'..') == 1
            mice_list_virus(o) = [];
        end
    end
end
Nmice_virus = length(mice_list_virus);
              

%% Import individual workspaces and load the data into a pooled array
%AnimalIDs
for nummice=1:length(mice_list_virus)
    AnimalIDs(nummice) = {['ID',mice_list_virus(nummice).name(end-2:end)]};
end  
PooledAnimalID_only = AnimalIDs;
%SessionIDs
for s=1:numberofsessions
    SessionIDs(s) = {['SessionNum',num2str(s)]};
end
% ANIMAL ID
for nummice=1:length(mice_list_virus)
    for s=1:numberofsessions
        PooledAnimalID.(AnimalIDs{nummice}).(SessionIDs{s}) = {};
    end
end

%% Loop
for nummice=1:length(mice_list_virus)
    % Define the path to the mouse and find the folders to analyze: 
    path2mouse = [PATH2DATA,'\',mice_list_virus(nummice).name,'\']; 
    % if there are several sessions inside the mouse folder
    % if you are looking for workspaces .mat
%     sessions = dir(fullfile(path2mouse,'*.mat'));

    % if you are looking for folder
    sessions = dir(path2mouse);
        %remove non relevant folders
        for o = length(sessions):-1:1
            if sessions(o).isdir == 0  %remove non-folders
                sessions(o) = [];
            elseif strcmp(sessions(o).name,'data') == 1 || strcmp(sessions(o).name,'figures') == 1 || strcmp(sessions(o).name,'results') == 1 ...
                 || strcmp(sessions(o).name,'.') == 1 || strcmp(sessions(o).name,'..') == 1 || contains(sessions(o).name,'.csv') == 1
                 sessions(o) = [];
            end
        end
        if isempty(sessions)
            sessions(1).name = []; %to have an existent session
        end

    %SessionIDs
    for s=1:length(sessions)
        SessionNames(s) = {sessions(s).name};   
    end       

    % Create a folder to save the data for this mouse and define the path
    if exist([PATH2SAVEFOLDER,'\',mice_list_virus(nummice).name,'\'],'dir') == 0
        mkdir([PATH2SAVEFOLDER,'\',mice_list_virus(nummice).name,'\']);
    end
    path2save_mouse = [PATH2SAVEFOLDER,'\',mice_list_virus(nummice).name,'\'];

    % Loop for all the sessions for the mouse: Sessions are experimental days you want to average together (replicates of each other) .........
    for s = 1:length(sessions) 
  

          done=0;
        % load the mouse/session matlab space
        load([path2mouse,sessions(s).name,'\IndividualData.mat']);

        % Initialization of pooled structure with individual trials   
        if nummice == 1 && s == 1
            % INDIVIDUAL TRIALS OF EACH MICE FOR ALL MICE
            PooledINDIV.(AnimalIDs{nummice}).(SessionIDs{s}).dFF_GPe = ones(1,length(time_vect))*nan;
            PooledINDIV.(AnimalIDs{nummice}).(SessionIDs{s}).dFF_SNr = ones(1,length(time_vect))*nan;
            PooledINDIV.(AnimalIDs{nummice}).(SessionIDs{s}).time_vect = ones(1,length(time_vect))*nan;  
            PooledINDIV.(AnimalIDs{nummice}).(SessionIDs{s}).TTL = ones(1,length(time_vect))*nan;  
            Allmice_RotarodOn = ones(length(mice_list_virus),numberofsessions)*nan;  
            Allmice_RotarodOff = ones(length(mice_list_virus),numberofsessions)*nan;  
            Pearson_Lag.(AnimalIDs{nummice}).pre = ones(numberofsessions,2039)*nan;
            Pearson_Lag.(AnimalIDs{nummice}).during = ones(numberofsessions,2039)*nan;
            Pearson_Lag.(AnimalIDs{nummice}).post = ones(numberofsessions,2039)*nan;
            
            %AUC
            for ku = 1:length(dFF_names)
                AUCtrials.(dFF_names{ku}) = ones(length(SessionIDs),length(mice_list_virus)*3)*nan; %3 for pre/during/post
                AUCaverage.(dFF_names{ku}) = ones(length(mice_list_virus),3)*nan; %5 mice, 3 for pre/during/post

                %MEAN BASELINE
                Baselinetrials.(dFF_names{ku}) = ones(length(SessionIDs),length(mice_list_virus)*3)*nan; %3 for pre/during/post
                Baselineaverage.(dFF_names{ku}) = ones(length(mice_list_virus),3)*nan; %5 mice, 3 for pre/during/post            
            end
            
            %PEARSON R
            PEARSONtrials = ones(length(SessionIDs),length(mice_list_virus)*3)*nan; %3 for pre/during/post
            PEARSONaverage = ones(length(mice_list_virus),3)*nan; %5 mice, 3 for pre/during/post
            
        end

        % Add data one animal at a time
        % ANIMAL ID and SESSION ID
        PooledAnimalID.(AnimalIDs{nummice}).(SessionIDs{s}) = SessionNames(s);
        
        % STORE INDIVIDUAL TRIALS
        PooledINDIV.(AnimalIDs{nummice}).(SessionIDs{s}).GPe = streams.ZScoredFF.GPe;
        PooledINDIV.(AnimalIDs{nummice}).(SessionIDs{s}).SNr = streams.ZScoredFF.SNr;
        PooledINDIV.(AnimalIDs{nummice}).(SessionIDs{s}).time_vect = time_vect;  
        
        % STORE PEARSON LAG
        Pearson_Lag.(AnimalIDs{nummice}).pre(s,:) = pearson_lag2plot_pre;
        Pearson_Lag.(AnimalIDs{nummice}).during(s,:) = pearson_lag2plot_during;
        Pearson_Lag.(AnimalIDs{nummice}).post(s,:) = pearson_lag2plot_post;
        
        % TTL
        startrod_time = epocs.TTLstart_time; 
        endrod_time = epocs.TTLstop_time;
        rotarodevents = [startrod_time;endrod_time];
        tmp = find(rotarodevents(1) <time_vect);
        ix_rotarod1 = tmp(1); 
        tmp = find(rotarodevents(2) <time_vect);
        ix_rotarod2 = tmp(1);
        TTL = zeros(1,length(time_vect)); 
        TTL(ix_rotarod1:ix_rotarod2) = 5;
        PooledINDIV.(AnimalIDs{nummice}).(SessionIDs{s}).TTL = TTL;  
        Allmice_RotarodOn(nummice,s) = startrod_time(1);  
        Allmice_RotarodOff(nummice,s) = endrod_time(1);  
        
        for ku = 1:length(dFF_names)
            %AUC
            AUCtrials.(dFF_names{ku})(s,nummice) = pre_auc.(dFF_names{ku});
            AUCtrials.(dFF_names{ku})(s,nummice+length(mice_list_virus)) = during_auc.(dFF_names{ku});
            AUCtrials.(dFF_names{ku})(s,nummice+2*length(mice_list_virus)) = post_auc.(dFF_names{ku});

            %Baseline
            Baselinetrials.(dFF_names{ku})(s,nummice) = MeanBasalLevelpre.(dFF_names{ku});
            Baselinetrials.(dFF_names{ku})(s,nummice+length(mice_list_virus)) = MeanBasalLevelduring.(dFF_names{ku});
            Baselinetrials.(dFF_names{ku})(s,nummice+2*length(mice_list_virus)) = MeanBasalLevelpost.(dFF_names{ku});
        end
        
        %PEARSON R
        PEARSONtrials(s,nummice) = corr_pre; %3 for pre/during/post
        PEARSONtrials(s,nummice+length(mice_list_virus)) = corr_during; %3 for pre/during/post
        PEARSONtrials(s,nummice+2*length(mice_list_virus)) = corr_post; %3 for pre/during/post

    end
end
%% end    

for ku = 1:length(dFF_names)
    % AUCaverage
    for i=1:length(mice_list_virus)
        AUCaverage.(dFF_names{ku})(i,1) = mean(AUCtrials.(dFF_names{ku})(1:length(SessionIDs),i),1);                                %pre
        AUCaverage.(dFF_names{ku})(i,2) = mean(AUCtrials.(dFF_names{ku})(1:length(SessionIDs),i+length(mice_list_virus)),1);        %during
        AUCaverage.(dFF_names{ku})(i,3) = mean(AUCtrials.(dFF_names{ku})(1:length(SessionIDs),i+2*length(mice_list_virus)),1);      %post
    end

    % Baselineaverage
    for i=1:length(mice_list_virus)
        Baselineaverage.(dFF_names{ku})(i,1) = mean(Baselinetrials.(dFF_names{ku})(1:length(SessionIDs),i),1);                              %pre
        Baselineaverage.(dFF_names{ku})(i,2) = mean(Baselinetrials.(dFF_names{ku})(1:length(SessionIDs),i+length(mice_list_virus)),1);      %during
        Baselineaverage.(dFF_names{ku})(i,3) = mean(Baselinetrials.(dFF_names{ku})(1:length(SessionIDs),i+2*length(mice_list_virus)),1);    %post
    end
end

% Pearsonaverage
for i=1:length(mice_list_virus)
    PEARSONaverage(i,1) = mean(PEARSONtrials(1:length(SessionIDs),i),1);                              %pre
    PEARSONaverage(i,2) = mean(PEARSONtrials(1:length(SessionIDs),i+length(mice_list_virus)),1);      %during
    PEARSONaverage(i,3) = mean(PEARSONtrials(1:length(SessionIDs),i+2*length(mice_list_virus)),1);    %post
end


% Lag
for i=1:length(mice_list_virus)
    pearson_lag_pre(i,:) = mean(Pearson_Lag.(AnimalIDs{i}).pre,1);                              %pre
    pearson_lag_during(i,:) = mean(Pearson_Lag.(AnimalIDs{i}).during,1);      %during
    pearson_lag_post(i,:) = mean(Pearson_Lag.(AnimalIDs{i}).post,1);    %post
end


%% Create new structure with data aligned at the "rotarod on" and all same length (add NaNs)

for ku = 1:length(dFF_names)

    % determine maximum length of structure
    length_all = ones(length(SessionIDs),length(mice_list_virus))*nan;
    for nummice=1:length(mice_list_virus)
         for s = 1:length(sessions)
            length_all(s,nummice) = length(PooledINDIV.(AnimalIDs{nummice}).(SessionIDs{s}).time_vect);
         end
    end
    max_length = max(length_all,[],'all');

    %new time vect
    max_idx = max_length;
    totaltimesec = dt_ds*max_idx; %number of seconds in the entire data
    time_vect_new = 0:dt_ds*1:dt_ds*(max_idx-1);

    % initialize new structure and fill it with nans
    for nummice=1:length(mice_list_virus)
        PooledINDIV_aligned.(dFF_names{ku}).(AnimalIDs{nummice}) = ones(length(SessionIDs),max_length)*nan;
    end

    % find the timepoint of rotarod start for all mice/sessions
    TTLon_all = ones(length(SessionIDs),length(mice_list_virus))*nan;
    TTLoff_all = ones(length(SessionIDs),length(mice_list_virus))*nan;

    for nummice=1:length(mice_list_virus)
         for s = 1:length(sessions)
            temp1 = find(PooledINDIV.(AnimalIDs{nummice}).(SessionIDs{s}).TTL);
            TTLon_all(s,nummice) = temp1(1);
            temp2 = find(PooledINDIV.(AnimalIDs{nummice}).(SessionIDs{s}).TTL);
            TTLoff_all(s,nummice) = temp2(length(temp2));
         end
    end
    earliestTTLon = min(TTLon_all,[],'all');

    % create for loop
    for nummice=1:length(mice_list_virus)
         for s = 1:length(sessions)
             length_dFF_temp = length(PooledINDIV.(AnimalIDs{nummice}).(SessionIDs{s}).time_vect);
             if TTLon_all(s,nummice) > earliestTTLon
                deltashift = TTLon_all(s,nummice) - earliestTTLon;
                PooledINDIV_aligned.(dFF_names{ku}).(AnimalIDs{nummice})(s,1:length_dFF_temp-deltashift) = PooledINDIV.(AnimalIDs{nummice}).(SessionIDs{s}).(dFF_names{ku})(1+deltashift:end);
                TTLon_all(s,nummice) = TTLon_all(s,nummice) - deltashift; %updated TTL on
                TTLoff_all(s,nummice) = TTLoff_all(s,nummice) - deltashift; % updated TTL off
             elseif TTLon_all(s,nummice) == earliestTTLon
                PooledINDIV_aligned.(dFF_names{ku}).(AnimalIDs{nummice})(s,1:length_dFF_temp) = PooledINDIV.(AnimalIDs{nummice}).(SessionIDs{s}).(dFF_names{ku});
             end
         end
    end

end

%% Create new structure with data aligned at the "rotarod OFF" and all same length (add NaNs)

for ku = 1:length(dFF_names)

    % initialize new structure and fill it with nans
    for nummice=1:length(mice_list_virus)
        PooledINDIV_alignedOFF.(dFF_names{ku}).(AnimalIDs{nummice}) = ones(length(SessionIDs),max_length)*nan;
    end

    % find the timepoint of rotarod STOP for all mice/sessions
    TTLon_all_4aligned2OFF = ones(length(SessionIDs),length(mice_list_virus))*nan;
    TTLoff_all_4aligned2OFF = ones(length(SessionIDs),length(mice_list_virus))*nan;

    for nummice=1:length(mice_list_virus)
         for s = 1:length(sessions)
            temp1 = find(PooledINDIV.(AnimalIDs{nummice}).(SessionIDs{s}).TTL);
            TTLon_all_4aligned2OFF(s,nummice) = temp1(1);
            temp2 = find(PooledINDIV.(AnimalIDs{nummice}).(SessionIDs{s}).TTL);
            TTLoff_all_4aligned2OFF(s,nummice) = temp2(length(temp2));
            clear temp1 temp2
         end
    end
    latestTTLoff = max(TTLoff_all_4aligned2OFF,[],'all');

    % create for loop
    for nummice=1:length(mice_list_virus)
         for s = 1:length(sessions)
             length_dFF_temp = length(PooledINDIV.(AnimalIDs{nummice}).(SessionIDs{s}).time_vect);
             if TTLoff_all_4aligned2OFF(s,nummice) < latestTTLoff
                deltashift2 = latestTTLoff - TTLoff_all_4aligned2OFF(s,nummice);
                if (length_dFF_temp + deltashift2) < max_length +1;              
                    PooledINDIV_alignedOFF.(dFF_names{ku}).(AnimalIDs{nummice})(s,deltashift2:length_dFF_temp+deltashift2-1) = PooledINDIV.(AnimalIDs{nummice}).(SessionIDs{s}).(dFF_names{ku});
                else
                    PooledINDIV_alignedOFF.(dFF_names{ku}).(AnimalIDs{nummice})(s,deltashift2:max_length) = PooledINDIV.(AnimalIDs{nummice}).(SessionIDs{s}).(dFF_names{ku})(1:max_length-deltashift2+1);                    
                end    
                TTLon_all_4aligned2OFF(s,nummice) = TTLon_all_4aligned2OFF(s,nummice) + deltashift2; %updated TTL on
                TTLoff_all_4aligned2OFF(s,nummice) = TTLoff_all_4aligned2OFF(s,nummice) + deltashift2; % updated TTL off
             elseif TTLoff_all_4aligned2OFF(s,nummice) == latestTTLoff
                PooledINDIV_alignedOFF.(dFF_names{ku}).(AnimalIDs{nummice})(s,1:length_dFF_temp) = PooledINDIV.(AnimalIDs{nummice}).(SessionIDs{s}).(dFF_names{ku});
             end
         end
    end

end






%% Save pooled data for all viruses together
% STORE INDIVIDUAL and AVERAGE TRIALS
save([PATH2SAVEFOLDER,'\PooledAllMice.mat'],'PooledINDIV_aligned','PooledINDIV_alignedOFF','AUCtrials','AUCaverage','PooledAnimalID','PooledAnimalID_only',...
    'TTLon_all_4aligned2OFF','TTLoff_all_4aligned2OFF','TTLoff_all','TTLon_all','Baselinetrials','Baselineaverage','PEARSONtrials','PEARSONaverage');



%% Plot the results per trial types separately
for ku=1:length(dFF_names)
    
avg_mode = 1; %1 if mean, 2 if median
color2plot = {'k','g','m','c','b','r'};


figure;        
if show_plot == 0
   set(gcf,'visible','off')
else
   set(gcf,'visible','on')
end

totalnum = length(mice_list_virus)*size(PooledINDIV_aligned.(dFF_names{ku}).(AnimalIDs{nummice}),1);
totallist = 2:1:totalnum+1;
c=gray(totalnum); %for more options see colormap; eg: spring summer jet cool hot parula hsv autumn winter gray bone copper pink lines

% optional: adjust time vect so that rotarod on is time 0
time_vect_adj = time_vect_new-round(min(TTLon_all,[],'all')*dt_ds); 

i=0;
for nummice=1:length(mice_list_virus)
    for s=1:size(PooledINDIV_aligned.(dFF_names{ku}).(AnimalIDs{nummice}),1)
%         i=round(rand*totalnum);
        i=i+1;
        plot(time_vect_adj,PooledINDIV_aligned.(dFF_names{ku}).(AnimalIDs{nummice})(s,:),'Color',[c(i,1) c(i,2) c(i,3)]); 
        hold on % plot individual trials
    end
end


if avg_mode == 1
    plot(time_vect_adj,nanmean(PooledINDIV_aligned.(dFF_names{ku}).(AnimalIDs{nummice}),1),'k',...
        'LineWidth',2) % plot the average on top of the individual trials
elseif avg_mode == 2
    plot(time_vect_adj,nanmedian(PooledINDIV_aligned.(dFF_names{ku}).(AnimalIDs{nummice}),1),'k',...
        'LineWidth',2)
end

%plot the rotarod area
rotarodarea = round((mean(TTLoff_all,'all') - mean(TTLon_all,'all'))*dt_ds)
error_area_onlyrectangle([0 (mean(TTLoff_all_4aligned2OFF,'all')*dt_ds-180)],[-4 8],[12 12],[0.8500, 0.1250, 0.0980],0.1); %t_trials is the time vector 

ax = gca;
ax.FontSize = 20; 
set(ax, 'box','off');
yline(0,'-.k');
xlabel('Time (s)','FontSize', 20);
ylabel('Zscore dFF (All Trials)','FontSize', 20);
title(['Rotarod On - ',dFF_names{ku}],'FontSize', 16);  
min_x=-35; max_x=300; min_y=-4; max_y=8; axis([min_x max_x min_y max_y]); 
set(gca,'TickLength',[0 0])
xticks([-30 0 30 60 90 120 150 180 210 240 270 300]);

% sgtitle([raw_or_corr{r},' ',datatype{i},' ',dFF_names{d}]); % for groups of subplots
if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
    saveas(gcf,[PATH2SAVEFOLDER,'\pooled figures\ALL INDIV TRIALS ',dFF_names{ku},'.tif'])
    saveas(gcf,[PATH2SAVEFOLDER,'\pooled figures\ALL INDIV TRIALS ',dFF_names{ku},'.fig'])
end

end


%% Plot the results per trial types separately, aligned to off



for ku=1:length(dFF_names)   
    avg_mode = 1; %1 if mean, 2 if median
    color2plot = {'k','g','m','c','b','r'};
    figure;        
    if show_plot == 0
       set(gcf,'visible','off')
    else
       set(gcf,'visible','on')
    end

    totalnum = length(mice_list_virus)*size(PooledINDIV_alignedOFF.(dFF_names{ku}).(AnimalIDs{nummice}),1);
    totallist = 2:1:totalnum+1;
    c=gray(totalnum); %

    % time
    TTLoff_all_4aligned2OFF_sec = mean(TTLoff_all_4aligned2OFF,'all')*dt_ds;
    TTLon_all_4aligned2OFF_sec = min(TTLon_all_4aligned2OFF,[],'all')*dt_ds;
    xt = round(TTLoff_all_4aligned2OFF_sec - TTLon_all_4aligned2OFF_sec);
  
    % optional: adjust time vect so that rotarod on is time 0
    time_vect_adj2 = time_vect_new-TTLon_all_4aligned2OFF_sec; 

    i=0;
    for nummice=1:length(mice_list_virus)
        for s=1:size(PooledINDIV_alignedOFF.(dFF_names{ku}).(AnimalIDs{nummice}),1)
    %         i=round(rand*totalnum);
            i=i+1;
            plot(time_vect_adj2,PooledINDIV_alignedOFF.(dFF_names{ku}).(AnimalIDs{nummice})(s,:),'Color',[c(i,1) c(i,2) c(i,3)]); 
            hold on % plot individual trials
        end
    end

    if avg_mode == 1
        plot(time_vect_adj2,nanmean(PooledINDIV_alignedOFF.(dFF_names{ku}).(AnimalIDs{nummice}),1),'k',...
            'LineWidth',2) % plot the average on top of the individual trials
    elseif avg_mode == 2
        plot(time_vect_adj2,nanmedian(PooledINDIV_alignedOFF.(dFF_names{ku}).(AnimalIDs{nummice}),1),'k',...
            'LineWidth',2)
    end
     %plot the rotarod area
     error_area_onlyrectangle([(TTLon_all_4aligned2OFF_sec-TTLon_all_4aligned2OFF_sec) (TTLoff_all_4aligned2OFF_sec-TTLon_all_4aligned2OFF_sec)],[-10 -10],[xt xt],[0.8500, 0.1250, 0.0980],0.1); %t_trials is the time vector 

    xline(TTLoff_all_4aligned2OFF_sec-TTLon_all_4aligned2OFF_sec,':','lineWidth',2,'Color','b');
    xline(TTLon_all_4aligned2OFF_sec-TTLon_all_4aligned2OFF_sec,':','lineWidth',2,'Color','b');

    ax = gca;
    ax.FontSize = 20; 
    set(ax, 'box','off');
    yline(0,'-.k');
    xlabel('Time (s)','FontSize', 20);
    ylabel('Zscore dFF (All Trials)','FontSize', 20);
    title(['Rotarod On, aligned to Off - ',dFF_names{ku}],'FontSize', 16);  
    min_x=-125; max_x=425; min_y=-4; max_y=10; axis([min_x max_x min_y max_y]); 
    set(gca,'TickLength',[0 0])
    xticks([-120 0 120 240 360]);

    % sgtitle([raw_or_corr{r},' ',datatype{i},' ',dFF_names{d}]); % for groups of subplots
    if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
        saveas(gcf,[PATH2SAVEFOLDER,'\pooled figures\ALL INDIV TRIALS ALIGNED2OFF',dFF_names{ku},'.tif'])
        saveas(gcf,[PATH2SAVEFOLDER,'\pooled figures\ALL INDIV TRIALS ALIGNED2OFF',dFF_names{ku},'.fig'])
    end
end




%% Create plot with nanmean and nanstd, average per mouse used for full graph _ same but also align to rotarod off (align to both)

for ku=1:length(dFF_names)
    color2plot = {[0.68,0.03,0.89],[1.00,0.00,0.00]};
    %plot
    figure; clf; hold on
    if show_plot == 0
        set(gcf,'visible','off')
    else
       set(gcf,'visible','on')
    end

    % ON
    maxlength_on = max(length(PooledINDIV_aligned.(dFF_names{ku}).(AnimalIDs{nummice})));
    %plot aligned to rotarod ON
    averagedFF_pre.(dFF_names{ku}) = ones(length(mice_list_virus),maxlength_on)*nan;
    for nummice=1:length(mice_list_virus)
        averagedFF_pre.(dFF_names{ku})(nummice,:) = nanmean(PooledINDIV_aligned.(dFF_names{ku}).(AnimalIDs{nummice}),1);
    end
    midrod1 = round(mean(TTLon_all,'all')) + 60./dt_ds; % 
    averagedFF_pre.(dFF_names{ku})(:,midrod1:end)=[];
    %step 1: normal timevect
    time_vect_pre = time_vect_new(1:midrod1); %
    %step 2: time vect starts at 0 when rotarod starts
    % time_vect_pre2 = time_vect_pre-round(min(TTLon_all,[],'all')*dt_ds); 
    time_vect_pre2 = time_vect_pre-(round(mean(TTLon_all,'all'))*dt_ds); 

 

    average_pool_pre.(dFF_names{ku}) = nanmean(averagedFF_pre.(dFF_names{ku}),1); %or slice even more
    error_pool_pre.(dFF_names{ku}) = nanstd(averagedFF_pre.(dFF_names{ku}),1,1)./sqrt(sum(~isnan(averagedFF_pre.(dFF_names{ku})(:,1))));
    error_area(time_vect_pre2,average_pool_pre.(dFF_names{ku}),error_pool_pre.(dFF_names{ku}),color2plot{ku},0.3); %t_trials is the time vector 

    % OFF
    length_off = ones(1,7);
    for nummice=1:length(mice_list_virus)
        length_off(nummice)=length(PooledINDIV_alignedOFF.(dFF_names{ku}).(AnimalIDs{nummice}));
    end
    maxlength_off=max(length_off);

    % % plot aligned to rotarod OFF
    averagedFF_post.(dFF_names{ku}) = ones(length(mice_list_virus),maxlength_off)*nan;
    for nummice=1:length(mice_list_virus)
        averagedFF_post.(dFF_names{ku})(nummice,:) = nanmean(PooledINDIV_alignedOFF.(dFF_names{ku}).(AnimalIDs{nummice}),1);
    end
    midrod2 = round(mean(TTLoff_all_4aligned2OFF,'all')) - 60./dt_ds;
    averagedFF_post.(dFF_names{ku})(:,1:midrod2)=[]; %or slice even more
    %step 1: normal timevect
    time_vect_post = time_vect_new; 
    time_vect_post(:,1:midrod2)=[]; 
    %step 2: time vect stops at 60 when rotarod stops
    time_vect_post2 = time_vect_post-round(min(TTLoff_all_4aligned2OFF,[],'all')*dt_ds)+120; %or slice even more
     %fill
    average_pool_post.(dFF_names{ku}) = nanmean(averagedFF_post.(dFF_names{ku}),1);
    error_pool_post.(dFF_names{ku})= nanstd(averagedFF_post.(dFF_names{ku}),1,1)./sqrt(sum(~isnan(averagedFF_post.(dFF_names{ku})(:,1))));
    error_area(time_vect_post2,average_pool_post.(dFF_names{ku}),error_pool_post.(dFF_names{ku}),color2plot{ku},0.3); %t_trials is the time vector 


      error_area_onlyrectangle([0 120],[-2 3],[120 120],[0.8500, 0.1250, 0.0980],0.1); %t_trials is the time vector 

      plot([60 60],[1.3 2.8],':','lineWidth',2,'Color','k');
    ax = gca;
    ax.FontSize = 20; 
    set(ax, 'box','off');
    yline(0,'-.k');
    xlabel('Time (s)','FontSize', 20);
    ylabel('Average Zscore dFF','FontSize', 20);
    title(['Rotarod On - ',dFF_names{ku}],'FontSize', 16);  
    min_x=-125; max_x=245; min_y=-4; max_y=5; axis([min_x max_x min_y max_y]); 
    set(gca,'TickLength',[0 0])
    xticks([-120 0 120 240]);

    %saveplot or not
    if done == 0 && save_plot == 1 
        saveas(gcf,[PATH2SAVEFOLDER,'\pooled figures\AVE all trials_doublealign',dFF_names{ku},'.tif']);
        saveas(gcf,[PATH2SAVEFOLDER,'\pooled figures\AVE all trials_doublealign',dFF_names{ku},'.fig']);
    end
end  


%% Create representative plot of 1 animal with GPe and SNr on the same plot

whichanimal = 3;
whichtrace = 4;

for ku=1:length(dFF_names)

% ON
%plot aligned to rotarod ON
selecteddFF_pre.(dFF_names{ku}) = PooledINDIV_aligned.(dFF_names{ku}).(AnimalIDs{whichanimal})(whichtrace,:);
midrod1 = round(mean(TTLon_all,'all')) + 60./dt_ds; % 
selecteddFF_pre.(dFF_names{ku})(:,midrod1:end)=[];
%step 1: normal timevect
time_vect_pre = time_vect_new(1:midrod1); %or slice even more
%step 2: time vect starts at 0 when rotarod starts
% time_vect_pre2 = time_vect_pre-round(min(TTLon_all,[],'all')*dt_ds); 
time_vect_pre2 = time_vect_pre-(round(mean(TTLon_all,'all'))*dt_ds); 

% OFF
% plot aligned to rotarod OFF
selecteddFF_post.(dFF_names{ku}) = PooledINDIV_alignedOFF.(dFF_names{ku}).(AnimalIDs{whichanimal})(whichtrace,:);
midrod2 = round(mean(TTLoff_all_4aligned2OFF,'all')) - 60./dt_ds;
selecteddFF_post.(dFF_names{ku})(:,1:midrod2)=[]; %or slice even more
%step 1: normal timevect
time_vect_post = time_vect_new; 
time_vect_post(:,1:midrod2)=[]; 
%step 2: time vect stops at 60 when rotarod stops
time_vect_post2 = time_vect_post-round(min(TTLoff_all_4aligned2OFF,[],'all')*dt_ds)+120; %or slice even more
end

%plot
figure; clf; hold on
if show_plot == 0
    set(gcf,'visible','off')
else
   set(gcf,'visible','on')
end

color2plot = {[1.00,0.00,0.00],[0.68,0.03,0.89]};

% plot
plot(time_vect_pre2,selecteddFF_pre.(dFF_names{1}),'Color',color2plot{1}); %t_trials is the time vector 
hold on
plot(time_vect_pre2,selecteddFF_pre.(dFF_names{2}),'Color',color2plot{2}); %t_trials is the time vector 
plot(time_vect_post2,selecteddFF_post.(dFF_names{1}),'Color',color2plot{1}); %t_trials is the time vector
plot(time_vect_post2,selecteddFF_post.(dFF_names{2}),'Color',color2plot{2}); %t_trials is the time vector

% plot properties
error_area_onlyrectangle([0 120],[-2 3],[120 120],[0.8500, 0.1250, 0.0980],0.1); %t_trials is the time vector 

% xline(mean(TTLoff_all_4aligned2OFF,'all')*dt_ds,':','lineWidth',2,'Color','b');
plot([60 60],[1.3 2.8],':','lineWidth',2,'Color','k');
ax = gca;
ax.FontSize = 20; 
set(ax, 'box','off');
yline(0,'-.k');
xlabel('Time (s)','FontSize', 20);
ylabel('Average Zscore dFF','FontSize', 20);
title(['Rotarod On'],'FontSize', 16);  
min_x=-125; max_x=245; min_y=-4; max_y=7; axis([min_x max_x min_y max_y]); 
set(gca,'TickLength',[0 0])
xticks([-120 0 120 240]);

%saveplot or not
if done == 0 && save_plot == 1 
    saveas(gcf,[PATH2SAVEFOLDER,'\pooled figures\Representative Trace ID 28,t4.tif']);
    saveas(gcf,[PATH2SAVEFOLDER,'\pooled figures\Representative Trace ID 28,t4.fig']);
end




%% Individual trials- heatmaps _ same but also align to rotarod off
for ku=1:length(dFF_names)

    
orderedmap = 0; % 
Color_scale = [0 3];
figure;clf
set(gcf,'visible','on')

if show_plot == 0
    set(gcf,'visible','off')
end

%prepare the plot
alltrials_array = [];
for nummice=1:length(mice_list_virus)
    animals = PooledAnimalID_only;
    % add data aligned to rotarod ON
    indivdata_pre = PooledINDIV_aligned.(dFF_names{ku}).(AnimalIDs{nummice}); %version 1: the whole thing
    % shift the time vector so that 0 is at rotarod start
    time_vect_pre = time_vect_new-round(mean(TTLon_all,'all')*dt_ds); 
    % version 2: crop anything past 60 sec
    index_crop=find((time_vect_pre-60)>0);
    index1 = index_crop(1);
    time_vect_pre(index1+1:end)=[];
    indivdata_pre(:,index1+1:end)=[];
    % add data aligned to rotarod OFF
    indivdata_post = PooledINDIV_alignedOFF.(dFF_names{ku}).(AnimalIDs{nummice}); %version 1: the whole thing
    % shift the time vector so that 0 is at rotarod start
    time_vect_post = time_vect_new-round(min(TTLoff_all_4aligned2OFF,[],'all')*dt_ds)+120; %or slice even more
    %crop: only start showing from timepoint 61
    index_crop2=find((time_vect_post-61)>0);
    index2 = index_crop2(1);
    time_vect_post(1:index2)=[];
    indivdata_post(:,1:index2)=[];
    nanbridge = ones(length(sessions),round(1./dt_ds))*nan;
    indivdata = [indivdata_pre nanbridge indivdata_post];
    time_vect_new2 = [time_vect_pre time_vect_post-1];
    numrows = size(indivdata,1);
    numcolumns = size(indivdata,2);
    if isempty(alltrials_array)                               
        alltrials_array(1:numrows,:) = indivdata;
    else
       numcurrentrows = size(alltrials_array,1);
       alltrials_array(numcurrentrows+1:numcurrentrows+numrows,:) = indivdata;
    end
end


% if you want each animal to have its trials sorted by order in terms of deg inhibition
if orderedmap == 1
    [~,order] = sort(abs(nanmean(alltrials_array(:,time_vect_new2 > 30 & time_vect_new2 < 90),2))); %if A is a matrix, then min(A,[],2) is a column vector containing the minimum value of each row
    %B = sort(A) sorts the elements of A in ascending order + [B,I] = sort(___) also returns a collection of index vectors
else
    order = [1:size(alltrials_array,1)]';
end
%plot
imagesc(time_vect_new2,1,alltrials_array(order,:))
colormap('copper')
%for more options see colormap; eg: spring summer jet cool hot parula hsv autumn winter gray bone copper pink lines


% Formatting
c = colorbar;
c.Label.String = 'Zscore dFF';
c.Label.FontSize = 20;
if ~isempty(Color_scale)
    c_limits = Color_scale;
    caxis(c_limits)
end
hold on
xlim([-125 245])
ax = gca;
ax.FontSize = 20; 
xlabel('Time (s)','FontSize', 20)
ylabel('Trials (All mice)','FontSize', 20)
tickarray = [0:10:size(alltrials_array,1)];
tickarray(1) = 1;
set(gca,'YTick',tickarray);
title(['Rotarod On - ',dFF_names{ku}],'FontSize', 16);
box off
ylimits = get(gca,'YLim');
plot([0 0],ylimits,'k','LineWidth',2,'LineStyle','-')
plot([120 120],ylimits,'k','LineWidth',2,'LineStyle','-')
plot([60 60],ylimits,'w','LineWidth',8,'LineStyle','-')
plot([61 61],ylimits,'k','LineWidth',2,'LineStyle',':')
plot([59 59],ylimits,'k','LineWidth',2,'LineStyle',':')
set(gca,'TickLength',[0 0])
xticks([-120 0 120 240]);

% clear alltrials_array
if save_plot == 1
    saveas(gcf,[PATH2SAVEFOLDER,'\pooled figures\Heatmaps of dFF signals_doublealign',dFF_names{ku},'.tif']);
    saveas(gcf,[PATH2SAVEFOLDER,'\pooled figures\Heatmaps of dFF signals_doublealign',dFF_names{ku},'.fig']);
end
                    
end               
