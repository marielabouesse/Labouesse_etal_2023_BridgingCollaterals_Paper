% Marie_FP_IndivData_extraction_GPeSNR_OF 
% DATA ANALYSIS PIPELINE FOR FIBER PHOTOMETRY EXPERIMENT: OPEN FIELD GPESNR: 
% Marie Labouesse 13 Jan 2022

% need to select the folder in which several animals data exist; in principle could have several sessions per mice but not implemented yes in detail
% extracts FP data (collected with TDT) and behavioral data (DLC) for individual mice and processes it
% aligns FP and behavioral data given that videoframes (used for DLC) were taken with a camera recorded in the TDT program (streams)
% corrects DLC data (upsample, interpolate), calculates mouse speed, generates open field map of mouse trajectory (new_pos)
% generates and saves correlation values and graphs between FP data and behavior for entire session (Ovrl_corr)
% identify trials triggered by onset and offset of movement (body_mov)
% generates and saves a matrix with individual trials aligned to specific event (Stim)
% generates and saves graphs aligned to specific event
% quantifies peak and latency to peak and minima and latency to minima for individual trials (Measurements, Measurements_con)
% the early parts of the code (dFF and DLC preprocessing) are mostly controlled in the first intro section of the code, but the later stuff has to be edited in the code itself

% functions needed:
% ovrl_corr_calculation_loop
% TDTbin2mat
% inpaint_nans

%% INITIALIZATIONS
close all; clear variables; clc;
set(0,'defaultfigurecolor',[1 1 1])

%% SETUP PARAMETERS
%%%%%%%%%% TTL EXTRACTION
ttl_extract = 'default'; %'default' or 'special' (needs to be 7 characters)
            
%%%%%%%%%% DATA EXTRACTION
channel_number = 2; % 1 if 1-site recording, 2 if 2-site recording
channel_names = {'SNr','GPe'}; % eg brain region
color_number = 2; % 2 if 405 and 465 recording, 3 if also 565
color_names = {'c405','c465'}; % can also be 565
dFF_number = 1; % how many dFF to analyze, can be 1 or 2, based on number of channels and colors (if more than 2, need to edit script)
dFF_names = {'SNr','GPe'}; %put them in this order: 1) channel 1, 465 and 2) channel 1, 560 OR channel 2: 465 (if other combinations, need to edit script)
% datatype = {'ZScoredFF','dFF','dFFsy','speed'}; %data analyses of interest -- remove speed if not analysing behavior
datatype = {'ZScoredFF','dFF'}; %data analyses of interest -- remove speed if not analysing behavior

%  BEHAVIOR
BehaviorDeeplabcut = 0;
BehaviorAnymaze = 0;
Behavior = 0;

%%%%%%%%%%% PREPROCESSING
% trimming
timetrim_start_set = 120; 
timetrim_end_set = 5; %
% low pass filtering of raw 405 and 470
lowpassfiltering = 1; % 1 if you want to low pass filter- typically yes- 0 for no
lowpassfreq = 3; % typically 1Hz
% detrending dFF
detrending = 1;% %% Write 1 for detrend dFF, normal detrend function of Matlab- typically no (Write 0) for open field
% substract 8th percentile as a baseline calculate using a moving window
deletemovingwindow = 1; % 1 to delete moving baseline 8th percentile, 0 if no.
% high pass filtering of dFF
highpassfiltering = 0; %'yes' if you want to high pass filter, 0 for no. Typically no for open field
highpassfreq = 0.005; % typically 0.005Hz

%%%%%%%%%%% 0; % 1 to loop, 0 if you only want to test one mouse and dont want to loop --> if so, edit the mouse number you want below: nummice_set
LoopOrNot = 1;
nummice_set = 5;
calccorrelations = 1;

%%%%%%%%%%% SHOW, SAVE, OVERWRITE PARAMETERS
show_plot = 1; % If 0, plots are not displayed
save_plot = 1; % If 0, plots are not saved
reanalysis = 1; % If 1, the code runs for sessions already analyzed. If 0, session excluded if analysed. This is based on the existence of the matlab space "IndividualData.mat"
overwrite = 1; % If reanalyzing data, 1 if you want to save them and overwrite results and 0 if not          
Collapse_sessions = 0; % If 0, analysis only of individual sessions. If 1, analysis of the group data (pooled)




%% DEFINE PATHS
path2data = uigetdir('select folder'); % SELECT FOLDER OF GROUP TO ANALYZE
mice_list = dir(path2data); %all things in this folder
% Define paths and data to analyze
path2savefolder = path2data; %Path to save


%% IDENTIFY MICE TO ANALYZE 
for o = length(mice_list):-1:1
    if mice_list(o).isdir == 0  %remove non-folders
        mice_list(o) = [];
    else
        if  strcmp(mice_list(o).name,'data') == 1 || strcmp(mice_list(o).name,'figures') == 1 ...   %remove folders with data or figures or ., ..
            || contains(mice_list(o).name,'data') || contains(mice_list(o).name,'figures') ...
            || contains(mice_list(o).name,'results') || contains(mice_list(o).name,'other') ...
            || strcmp(mice_list(o).name,'.') == 1 || strcmp(mice_list(o).name,'..') == 1
            mice_list(o) = [];
        end
    end
end
Nmice = length(mice_list);  %number of mice to analyze


%% LOOP ACROSS MICE
if LoopOrNot == 0  %if you only want to analyze one animal
    Nmice = 1;
end
for nummice = 1:Nmice   
    if LoopOrNot == 0
        nummice = nummice_set;  %need to set the number of the mouse you want to analyze up above in PARAMETERS
    end
    AnimalID = ['ID_',mice_list(nummice).name];   %last 4 digits of folders should be the animal name (end-3:end)
    Dirsplit = strsplit(path2data,'\');
    Virus_cell = Dirsplit(length(Dirsplit));
    Virus = Virus_cell{:}; %just extracting the string out of the cell array   
    BehavData = {};

    %% Setting the paths for the mouse and the sessions
    % Define the path to the mouse and find the folders to analyze:
    path2mouse = [path2data,'\',mice_list(nummice).name,'\'];
    % if there are several sessions inside the mouse folder
    sessions = dir(path2mouse);
    %remove non relevant folders
    for o = length(sessions):-1:1
        if sessions(o).isdir == 0  %remove non-folders
            sessions(o) = [];
        elseif strcmp(sessions(o).name,'data') == 1 || strcmp(sessions(o).name,'figures') == 1 || strcmp(sessions(o).name,'results') == 1 || strcmp(sessions(o).name,'anymaze') == 1 ...
             || strcmp(sessions(o).name,'.') == 1 || strcmp(sessions(o).name,'..') == 1 || contains(sessions(o).name,'.csv') == 1 ...
            || contains(sessions(o).name,'representative figures') ||  contains(sessions(o).name,'figures_paper_0195') ...
            sessions(o) = [];
        end
    end
    if isempty(sessions)
        sessions(1).name = []; %to have an existent session
    end
    
    % Create a folder to save the data for this mouse and define the path
    if exist([path2savefolder,'\',mice_list(nummice).name,'\'],'dir') == 0
        mkdir([path2savefolder,'\',mice_list(nummice).name,'\'])
    end
    path2save_mouse = [path2savefolder,'\',mice_list(nummice).name,'\'];
        
    
    %% Loop for all the sessions for the mouse
    for s = 1:length(sessions)
        % Define the path to the session and create folder to save if needed:
        PATH2SESSION = [path2mouse,sessions(s).name];
        if exist([path2save_mouse,sessions(s).name],'dir') == 0
            mkdir([path2save_mouse,sessions(s).name])
        end
        if  length(sessions) == 1
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
        done = exist([PATH2SAVE,'IndividualData',AnimalID,'.mat'],'file'); %if you specify 'file' then matlab searches for both files and folders
        if done == 0 || reanalysis == 1 %if not done or if wish to reanalyze

    

            %% Select FP data file and extract FP data streams 
            files_tev = dir(fullfile(PATH2SESSION,'*.tev'));
            trace_name = files_tev.name(1:end-4); %remove .mat or.tev
            data = TDTbin2mat(PATH2SESSION, 'TYPE', {'epocs', 'scalars', 'streams'});

            % Create 'streams' structure, which will contain extracted raw data from each channel/color
            for channel=1:length(channel_names)
                for colore=1:length(color_names)
                    streams.rawdata.(channel_names{channel}).(color_names{colore}) = [];
                end
            end
            % first channel
            if contains(PATH2SESSION,'inverted')
                stream_name_type405A = {'x05B','05B','B05B','x05C'};
                stream_name_type465A = {'x65B','65B','B65B','x65C'};
                stream_name_type560A = {'x60A','60A','A60A','x60C'};
            else
                stream_name_type405A = {'x05A','05A','A05A'}; 
                stream_name_type465A = {'x65A','65A','A65A'}; 
                stream_name_type560A = {'x60A','60A','A60A'};                
            end
            
            for i=1:length(stream_name_type405A)
                if isfield(data.streams, stream_name_type405A{i}); %isfield checks if 'c' is a field inside structure a.b.c: isfield(a,'c')
    %                     fieldnames(data.streams) %to get directly the field names in data.streams
                    streams.rawdata.(channel_names{1}).c405 = data.streams.(stream_name_type405A{i}).data; %access data in data.stream using the name of the 405 store
                    streams.rawdata.(channel_names{1}).c465 = data.streams.(stream_name_type465A{i}).data; %access data in data.stream using the name of the 465 store
                    % if also recording in red
                    if color_number == 3           
                        streams.rawdata.(channel_names{1}).c560 = data.streams.(stream_name_type560A{i}).data; %access data in data.stream using the name of the 560 store
                    end
                end
            end
            % second channel (if there is)
            if channel_number == 2
                if contains(PATH2SESSION,'inverted')
                    stream_name_type405B = {'x05A','05A','A05A'}; 
                    stream_name_type465B = {'x65A','65A','A65A'}; 
                    stream_name_type560B = {'x60A','60A','A60A'};
                else
                    stream_name_type405B = {'x05B','05B','B05B','x05C'};
                    stream_name_type465B = {'x65B','65B','B65B','x65C'};
                    stream_name_type560B = {'x60A','60A','A60A','x60C'};
                end
                
                for i=1:length(stream_name_type405B)
                    if isfield(data.streams, stream_name_type405B{i}); %isfield checks if 'c' is a field inside structure a.b.c: isfield(a,'c')
                        streams.rawdata.(channel_names{2}).c405 = data.streams.(stream_name_type405B{i}).data; %access data in data.stream using the name of the 405 store
                        streams.rawdata.(channel_names{2}).c465 = data.streams.(stream_name_type465B{i}).data; %access data in data.stream using the name of the 465 store
                        % if also recording in red
                        if color_number == 3           
                            streams.rawdata.(channel_names{2}).c560 = data.streams.(stream_name_type560B{i}).data; %access data in data.stream using the name of the 560 store
                        end
                    end
                end
            end
            %Calculate min length
            length_stream = [];
            for channel=1:length(channel_names)
                for colore=1:length(color_names)
                    length_stream(end+1) = length(streams.rawdata.(channel_names{channel}).(color_names{colore}));
                end
            end
            min_length_stream = min(length_stream);
            %Adjust the lengths of the channels
            for channel=1:length(channel_names)
                for colore=1:length(color_names)
                    if length(streams.rawdata.(channel_names{channel}).(color_names{colore})) > min_length_stream
                       streams.rawdata.(channel_names{channel}).(color_names{colore})...
                        (min_length_stream+1:length(streams.rawdata.(channel_names{channel}).(color_names{colore}))) = []; 
                    end
                end
            end        

            %% Downsampling FP data, create a new structure for it (streams) and define time vector
            % sampling rate
            if isfield (data.streams, 'A65A')
                sampling_rate = data.streams.A65A.fs;
            elseif isfield (data.streams, 'x65A' )
                sampling_rate = data.streams.x65A.fs;               
            end

            N = 10; %downsample 10 times 
            sampling_rate_ds = sampling_rate/N;
            % downsample
            for channel=1:length(channel_names)
                for colore=1:length(color_names)
                    streams.rawdata.(channel_names{channel}).(color_names{colore}) = downsample(streams.rawdata.(channel_names{channel}).(color_names{colore}),10);
                end
            end
            length_data = length(streams.rawdata.(channel_names{channel}).(color_names{colore})); 
            % time vector
            max_idx = length_data;
            dt_ds = 1/sampling_rate_ds; % delta t or unit of time (in seconds) per sample --> for 101.73 Hz sampling rate (after downsampling), should be close to 10 ms (9.8 ms)
            time_vect = 0:dt_ds*1:dt_ds*(max_idx-1); %time vector; values in seconds, same number of data points as the data


            %% Extract TTL data (epocs) for opto stim and camera frames
            % Initialize epocs structure
            if ttl_extract == 'default' 
                epoc_type = {'PulseStart_time','PulseStop_time','StimAmp','StimTime','TTLstart_time','TTLstop_time','Testsync_start_time','Testsync_stop_time','Framesync_start_time',...
                'Possync_start_time','FrontCam_start_time','BackCam_start_time'};
                epocs = struct();
                for i=1:length(epoc_type)
                    epocs.(epoc_type{i}) = [];
                end
                % Overall opto pulse
                if isfield(data.epocs,'Pu1_');
                    epocs.PulseStart_time = data.epocs.Pu1_.onset;
                    epocs.PulseStop_time = data.epocs.Pu1_.offset;
                end
                % Opto Stim amplitude
                if isfield(data.epocs,'Am1_');
                    epocs.StimAmp = data.epocs.Am1_.data;
                    epocs.StimTime = data.epocs.Am1_.onset;
                end
                % Generic TTL in PCO
                if isfield(data.epocs,'PC0_');
                    epocs.TTLstart_time = data.epocs.PC0_.onset;   
                    epocs.TTLstop_time = data.epocs.PC0_.offset;   
                end
                % Test start/stop sync TTLs arriving from Anymaze
                if isfield(data.epocs,'PC1_');
                    epocs.Testsync_start_time = data.epocs.PC1_.onset;
                    epocs.Testsync_stop_time = data.epocs.PC1_.offset;
                end
                
                 % Test start/stop sync TTLs arriving from Anymaze
                if isfield(data.epocs,'U11_');
                    epocs.TTLstart_time = data.epocs.U11_.onset;
%                     epocs.Testsync_stop_time = data.epocs.PC1_.offset;
                end
                
                 % Test start/stop sync TTLs arriving from Anymaze
                if isfield(data.epocs,'U21_');
                    epocs.TTLstop_time = data.epocs.U21_.onset;
%                     epocs.Testsync_stop_time = data.epocs.PC1_.offset;
                end
                
                % Cam Frame sync TTLs arriving from Anymaze
                if isfield(data.epocs,'PC2_');
                    epocs.Framesync_start_time = data.epocs.PC2_.onset;    %PC2 in new program
                end 
                % % Position sync TTLs arriving from Anymaze
                if isfield(data.epocs,'PC3_');
                    epocs.Possync_start_time = data.epocs.PC3_.onset;
                end

                % % Cam1 camera frames taken directly in Synapse
                if isfield(data.epocs,'Cam1');
                    epocs.FrontCam_start_time = data.epocs.Cam1.onset;
                end

                % % Cam2 camera frames taken directly in Synapse
                if isfield(data.epocs,'Cam2');
                    epocs.BackCam_start_time = data.epocs.Cam2.onset;
                end

                % if the TTLs are organized differently:
            elseif ttl_extract == 'special'
                epoc_type = epoc_names;
                for i=1:length(epoc_type)
                    epocs.(epoc_type{i}) = data.epocs.(epoc_locations{i}).(epoc_locations2{i});
                end
            end     

            
            %% ROTAROD specific
            startrod_time = epocs.TTLstart_time; 
            endrod_time = epocs.TTLstop_time;
             
            if isempty(endrod_time)
                if contains(PATH2SESSION,'endrod275')
                    endrod_time = startrod_time + 275;
                    epocs.TTLstop_time = endrod_time;
                elseif contains(PATH2SESSION,'endrod270')
                    endrod_time = startrod_time + 270;
                    epocs.TTLstop_time = endrod_time;
                end
            end
            
            
            %% Clear TDT data structure (now everything we need is in "streams" or "epocs")
%             clear data

            %% TRIMMING           
            % Setting up trim indexes
                timetrim_start = timetrim_start_set*1; %in seconds: adjust this manually
                timetrim_end =  timetrim_end_set*1; %in seconds: adjust this manually --> set up at the beginning of the code
                dummie1 = 1:length(time_vect); %indexes, same number as there are samples in time_vect
                dummie2 = dummie1(time_vect > timetrim_start); %only keep the indexes starting from the new start
                idx_start = dummie2(1); %index for the new start
                dummie2 = dummie1(time_vect > time_vect(end) - timetrim_end); %only keep the indexes from the new end to the current end
                idx_end = dummie2(1); %index for the new end
                clear dummie1 dummie2 
                      
            % Trimming and readjusting everything
            for channel=1:length(channel_names)
                for colore=1:length(color_names)
                    streams.rawdata_trim.(channel_names{channel}).(color_names{colore}) = streams.rawdata.(channel_names{channel}).(color_names{colore})(idx_start:idx_end);
                end
            end        
            time_vect_trim = time_vect(idx_start:idx_end)-time_vect(idx_start);
            %pulses
            epocs_trim = struct();
            for i=1:length(epoc_type)
                epocs_trim.(epoc_type{i}) = epocs.(epoc_type{i}) - timetrim_start;
            end
            
            % Trimming: Plot
            for channel=1:length(channel_names)
                for colore=1:length(color_names)
                    figure; clf; 
                    if show_plot == 0
                       set(gcf,'visible','off')
                    end
                    %Plot after
                    subplot(2,1,1);
                    plot(time_vect,streams.rawdata.(channel_names{channel}).(color_names{colore})); hold on; 
                    min_x=-200; max_x=max(time_vect)+200; min_y=mean(streams.rawdata.(channel_names{channel}).(color_names{colore}))-100; 
                    max_y=mean(streams.rawdata.(channel_names{channel}).(color_names{colore}))+100; axis([min_x max_x min_y max_y]); 
                    xlabel('time(sec)'); ylabel('fluorescence'); title('Data before trimming'); 
                    xline(timetrim_start,':k','trim'); xline(time_vect(end)-timetrim_end,':k','trim'); 
                    %Plot after
                    subplot(2,1,2); plot(time_vect_trim,streams.rawdata_trim.(channel_names{channel}).(color_names{colore})); hold on; 
                    axis([min_x max_x min_y max_y]); %ie same as other subplot
                    xlabel('time(sec)'); ylabel('fluorescence'); title('Data after trimming'); 
                    annotation('textbox',[.25 0 .1 .17],'string','Ru ok with the trimming? if not readjust','FitBoxToText','on','EdgeColor','none');
                    sgtitle([channel_names{channel},' ',color_names{colore}]);
                    if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                        saveas(gcf,[PATH2SAVE,'figures\trimming ',channel_names{channel},' ',color_names{colore},'.tif'])
                        saveas(gcf,[PATH2SAVE,'figures\trimming ',channel_names{channel},' ',color_names{colore},'.fig'])
                    end
                end 
            end

            % Reassign data variables if happy with trimming 
            for channel=1:length(channel_names)
                for colore=1:length(color_names)
                    streams.rawdata.(channel_names{channel}).(color_names{colore}) = streams.rawdata_trim.(channel_names{channel}).(color_names{colore});
                end
            end
            time_vect = time_vect_trim;
            for i=1:length(epoc_type)
                epocs.(epoc_type{i}) = epocs_trim.(epoc_type{i});
            end
           
            startrod_time = epocs.TTLstart_time; 
            endrod_time = epocs.TTLstop_time;
            
            length_data = length(streams.rawdata.(channel_names{channel}).(color_names{colore}));
            clear epocs_trim time_vect_trim 
            streams = rmfield(streams,'rawdata_trim'); %clear variable within structure.

                         

            %% LOW PASS FILTER OF FP DATA TO 1Hz
            if lowpassfiltering == 1
                %filter option 2: butter
                ftype = 'low';
                n_trials = 2; % 2nd order filter
                Wn = lowpassfreq/((sampling_rate_ds)/2); %lowpassfreq defined above
                % 0.5 Hz = 2 sec ; 1 Hz = 1 sec ; 2 Hz = 0.5 sec ; 3 Hz = 0.33 sec
                [a,b] = butter(n_trials,Wn,ftype);
                for channel=1:length(channel_names)
                    for colore=1:length(color_names)
                        streams.lowfilt.(channel_names{channel}).(color_names{colore}) = filtfilt(a,b,double(streams.rawdata.(channel_names{channel}).(color_names{colore})));
                        % plot
                        figure; clf; 
                        if show_plot == 0
                           set(gcf,'visible','off')
                        end
                        plot(time_vect,streams.rawdata.(channel_names{channel}).(color_names{colore})); hold on; 
                        plot(time_vect,streams.lowfilt.(channel_names{channel}).(color_names{colore}));
                        min_x=-10; max_x=max(time_vect)+10; min_y=min(streams.rawdata.(channel_names{channel}).(color_names{colore}))-1; 
                        max_y=max(streams.rawdata.(channel_names{channel}).(color_names{colore}))+1; axis([min_x max_x min_y max_y]); 
                        xlabel('samples'); ylabel('fluorescence'); title([channel_names{channel},' ', color_names{colore},' data before and after filtering']);
                        L=legend('data -10','data after filter'); L.Location = 'Best';
                        if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                            saveas(gcf,[PATH2SAVE,'figures\lowpass ',channel_names{channel},' ',color_names{colore},'.tif'])
                            saveas(gcf,[PATH2SAVE,'figures\lowpass ',channel_names{channel},' ',color_names{colore},'.fig'])
                        end
                    end
                end

                %% Save data if happy with low-pass filtering
                for channel=1:length(channel_names)
                    for colore=1:length(color_names)
                        streams.rawdata.(channel_names{channel}).(color_names{colore}) = streams.lowfilt.(channel_names{channel}).(color_names{colore});
                    end
                end
                streams = rmfield(streams,'lowfilt'); %clear variable within structure.
            end %of low pass filtering

                        
            

            %% CALCULATE dFF (synapse and polyfit versions)
            % Collect data from 'Pre' and 'Post' epochs to improve the polyfit (which allows to identify bleaching artefacts and correct for them)
            % rotarod ix
            rotarodevents = [startrod_time;endrod_time]
            tmp = find(rotarodevents(1) <time_vect);
            ix_rotarod1 = tmp(1); 
            tmp = find(rotarodevents(2) <time_vect);
            ix_rotarod2 = tmp(1); 
            
            % Calculations in each channel and color
            
            %% First channel, 465 color, first dFF
            channel=1;
            % dFF polyfit
            % Remove pre-post portions
            F465_1_prepost = [streams.rawdata.(channel_names{channel}).c465(1:ix_rotarod1),streams.rawdata.(channel_names{channel}).c465(ix_rotarod2:end)];
            F405_1_prepost = [streams.rawdata.(channel_names{channel}).c405(1:ix_rotarod1),streams.rawdata.(channel_names{channel}).c405(ix_rotarod2:end)];
            % Step 1: Polyfit on entire data 
            calc_coeff_fit_465 = polyfit(F405_1_prepost,F465_1_prepost,1); %fitting control Dv2 against signal Dv1 using least squares fit of degree 1
            % Step 2: Apply polyfit on 465 
            data405_fitted_465 = calc_coeff_fit_465(1).*streams.rawdata.(channel_names{channel}).c405 + calc_coeff_fit_465(2);
            % Step 3: Normalize
            streams.dFF.(dFF_names{1}) = 100*((streams.rawdata.(channel_names{channel}).c465 - data405_fitted_465)./data405_fitted_465); % of deltaF/F              

            % second dFF, either 2nd color in channel 1 or 2nd channel (1 color)- if ever I had 2 color, 2 channels, would need to edit
            % this script
            if color_number == 3 && channel_number == 1
                channel = 1;
                % dFF polyfit
                % Remove pre-post portions
                F560_1_prepost = [streams.rawdata.(channel_names{channel}).c560(1:ix_rotarod1),streams.rawdata.(channel_names{channel}).c560(ix_rotarod2:end)];
                % Step 1: Polyfit on entire data
                calc_coeff_fit_560 = polyfit(F405_1_prepost,F560_1_prepost,1); %fitting control Dv2 against signal Dv1 using least squares fit of degree 1
                % Step 2: Apply polyfit on 560 
                data405_fitted_560 = calc_coeff_fit_560(1).*streams.rawdata.(channel_names{channel}).c405 + calc_coeff_fit_560(2);
                % Step 3: Normalize
                streams.dFF.(dFF_names{2}) = 100*((streams.rawdata.(channel_names{channel}).c560 - data405_fitted_560)./data405_fitted_560); % of deltaF/F 
            elseif color_number == 2 && channel_number == 2
                channel = 2;
                % dFF polyfit
                % Remove pre-post portions
                F465_2_prepost = [streams.rawdata.(channel_names{channel}).c465(1:ix_rotarod1),streams.rawdata.(channel_names{channel}).c465(ix_rotarod2:end)];
                F405_2_prepost = [streams.rawdata.(channel_names{channel}).c405(1:ix_rotarod1),streams.rawdata.(channel_names{channel}).c405(ix_rotarod2:end)];
                % Step 1: Polyfit on entire data 
                calc_coeff_fit_465_2 = polyfit(F405_2_prepost,F465_2_prepost,1); %fitting control Dv2 against signal Dv1 using least squares fit of degree 1
                % Step 2: Apply polyfit on 465 
                data405_fitted_465_2 = calc_coeff_fit_465_2(1).*streams.rawdata.(channel_names{channel}).c405 + calc_coeff_fit_465_2(2);
                % Step 3: Normalize
                streams.dFF.(dFF_names{2}) = 100*((streams.rawdata.(channel_names{channel}).c465 - data405_fitted_465_2)./data405_fitted_465_2); % of deltaF/F  
            end

            % plot all
            for d=1:length(dFF_names)
                figure; clf; 
                if show_plot == 0
                   set(gcf,'visible','off')
                end
                yyaxis left
                plot(time_vect,streams.dFF.(dFF_names{d})); hold on; 
                ylabel('dFF');
                min_y=min(streams.dFF.(dFF_names{d}))-3; max_y=max(streams.dFF.(dFF_names{d}))+3; axis([min_x max_x min_y max_y]);
                yyaxis right
%                 plot(time_vect,streams.dFFsy.(dFF_names{d})); 
                for o = 1:size(epocs.PulseStart_time)
                    plot([epocs.PulseStart_time(o) epocs.PulseStart_time(o)],[-2 2],'m')
                    plot([epocs.PulseStop_time(o) epocs.PulseStop_time(o)],[-2 2],'c')
                end
                xlabel('time (sec)'); ylabel('dFF synapse'); title([dFF_names{d},' Normalized data']);
                L=legend('dFF','dFF Synapse','Stim on','Stim off'); L.Location = 'Best';
                min_x=-30; max_x=max(time_vect)+30; min_y=min(streams.dFF.(dFF_names{d}))-3; max_y=max(streams.dFF.(dFF_names{d}))+3; axis([min_x max_x min_y max_y]); 
                if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                    saveas(gcf,[PATH2SAVE,'figures\dFF calculations ',dFF_names{d},'.tif'])
                    saveas(gcf,[PATH2SAVE,'figures\dFF calculations ',dFF_names{d},'.fig'])
                end
            end

            %% Delete raw data (405 and 465) within streams
            streams = rmfield(streams,'rawdata'); %clear variable within structure.
            clear F0_405 F0_465 dFF_405 dFF_465 F0_405_2 F0_465_2 dFF_405_2 dFF_465_2 F0_560 dFF_560

            %% DETRENDING dFF
            if detrending == 1;% %% Detrend dFF, Option 1: normal detrend (other detrend options not here, but could include them later, perhaps as a function)
                for d=1:length(dFF_names)
                    streams.dFF_dtr.(dFF_names{d}) = detrend(streams.dFF.(dFF_names{d})); % 
                    %plot
                    figure; clf
                    if show_plot == 0
                       set(gcf,'visible','off')
                    end
                    h1 = plot(time_vect,streams.dFF.(dFF_names{d})); hold on; h2 = plot(time_vect,streams.dFF_dtr.(dFF_names{d}));
                    for o = 1:size(epocs.PulseStart_time)
                        plot([epocs.PulseStart_time(o) epocs.PulseStart_time(o)],[-2 2],'m')
                        plot([epocs.PulseStop_time(o) epocs.PulseStop_time(o)],[-2 2],'c')
                    end
                    min_x=-30; max_x=max(time_vect)+30; min_y=min(streams.dFF.(dFF_names{d}))-10; max_y=max(streams.dFF.(dFF_names{d}))+10; axis([min_x max_x min_y max_y]); 
                    xlabel('time (sec)'); ylabel('dFF'); title([dFF_names{d},' Normalized data, Regular detrend']);
                    L=legend([h1 h2],'dFF','dFF detrended'); L.Location='Best'; 
                    if done == 0 && save_plot == 1 || done == 0 && reanalysis == 1 && overwrite == 1
                        saveas(gcf,[PATH2SAVE,'figures\dFF detrending ',dFF_names{d},'.tif'])
                        saveas(gcf,[PATH2SAVE,'figures\dFF detrending ',dFF_names{d},'.fig'])
                    end
                    % if happy with detrend:
                    streams.dFF.(dFF_names{d}) = streams.dFF_dtr.(dFF_names{d});
                    streams = rmfield(streams,'dFF_dtr'); %clear variable within structure.
                end
            else
            end

           
            
            %% HIGH-PASS FILTERING
            if highpassfiltering == 1; %'yes' if you want to high pass filter- typically off for open field
                % Filter out the slow fluctuation
                ftype = 'high';
                n_trials = 2; % 2nd order filter
                Wn = highpassfreq/((sampling_rate_ds)/2); %highpassfreq: 
                [a,b] = butter(n_trials,Wn,ftype);

                for d=1:length(dFF_names)
                    streams.dFF_hp.(dFF_names{d}) = filtfilt(a,b,double(streams.dFF.(dFF_names{d})));
                    % plot
                    figure; clf
                    if show_plot == 0
                       set(gcf,'visible','off')
                    end
                    yyaxis left
                    h1=plot(time_vect,streams.dFF.(dFF_names{d})); hold on; 
                    for o = 1:size(epocs.PulseStart_time)
                        plot([epocs.PulseStart_time(o) epocs.PulseStart_time(o)],[-2 2],'m')
                        plot([epocs.PulseStop_time(o) epocs.PulseStop_time(o)],[-2 2],'c')
                    end
                    yyaxis right
                    h2=plot(time_vect,streams.dFF_hp.(dFF_names{d}));
                    min_x=-10; max_x=max(time_vect)+10; min_y=min(streams.dFF.(dFF_names{d}))-10; max_y=max(streams.dFF.(dFF_names{d}))+10; axis([min_x max_x min_y max_y]); 
                    xlabel('samples'); ylabel('fluorescence'); title([dFF_names{d},' dFF before and after filtering']);
                    L=legend([h1,h2],'dFF','dFF after high-pass filter'); L.Location = 'Best';

                    if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                        saveas(gcf,[PATH2SAVE,'figures\dFF highpass ',dFF_names{d},'.tif'])
                        saveas(gcf,[PATH2SAVE,'figures\dFF highpass ',dFF_names{d},'.fig'])
                    end
                    % if happy with detrend:
                    streams.dFF.(dFF_names{d}) = streams.dFF_hp.(dFF_names{d});
                    streams = rmfield(streams,'dFF_hp'); %clear variable within structure.
                end
            else
            end


            %% CALCULATE ROBUST Z-SCORES
            % Calculations
            % Calculate median, Median Absolute Deviation (MAD) and regular "modified" Z-score 
            for d=1:length(dFF_names)
                med_dFF.(dFF_names{d}) = median(streams.dFF.(dFF_names{d}));    
                MAD_dFF.(dFF_names{d}) = mad(streams.dFF.(dFF_names{d}),1);     
                streams.ZScoredFF.(dFF_names{d}) = 0.6745*(streams.dFF.(dFF_names{d})-med_dFF.(dFF_names{d}))./MAD_dFF.(dFF_names{d}); 
                % Plot
                figure; clf; 
                if show_plot == 0
                   set(gcf,'visible','off')
                end
                yyaxis left
                h1 = plot(time_vect,streams.dFF.(dFF_names{d}),'LineWidth',0.5,'Color','k'); hold on; 
                min_y = min(streams.dFF.(dFF_names{d}))-10; max_y = max(streams.dFF.(dFF_names{d}))+10; ylim([min_y max_y]);
                ylabel('dFF');
                yyaxis right
                h2 = plot(time_vect,streams.ZScoredFF.(dFF_names{d}), 'LineWidth',0.5,'Color',[1.0 0.2 0.2]);
                h3 = yline(0,'linewidth',2,'linestyle',':');
                for o = 1:size(epocs.PulseStart_time)
                    plot([epocs.PulseStart_time(o) epocs.PulseStart_time(o)],[-2 2],'m')
                    plot([epocs.PulseStop_time(o) epocs.PulseStop_time(o)],[-2 2],'c')
                end
                min_y = -10; max_y = 10; ylim([min_y max_y]);
                L = legend([h1,h2],'dFF','Z-score');
                L.Location = 'northeast';
                min_y = min(streams.ZScoredFF.(dFF_names{d}))-10; max_y = max(streams.ZScoredFF.(dFF_names{d}))+10;
                min_x = -20; max_x = max(time_vect)+20;  
                xlim([min_x max_x])
                xlabel('time (s)'); ylabel('Z-score dFF'); title([dFF_names{d},' dFF and Z-score']);
                if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                    saveas(gcf,[PATH2SAVE,'figures\dFF Zscore ',dFF_names{d},'.tif'])
                    saveas(gcf,[PATH2SAVE,'figures\dFF Zscore ',dFF_names{d},'.fig'])
                end
            end

          
             %% ADJUSTING dFF BASELINE TO 0 IN A MOVING WINDOW      
            if deletemovingwindow == 1;
                for d=1:length(dFF_names)
                    
                    % Define parameters for moving window 
                    baseline_window = 60; %ADJUST!
                    baseline_dt = 10; %Moving steps 
                    percentile = 8; %8th percentile
                    
                    % New time vector
                    temp_t = time_vect - time_vect(1); %creating new time vector set to starting at 0, incase time_vect weas not starting at 0
                    dummie = 1:length(temp_t); %generates index of same amount of samples (better than find cos too much computing time)
                    dummie = dummie(temp_t >= baseline_window); %value when the time is beyond the window --> so that we can apply the moving window calculations
                    idx_baseline_window = dummie(1); %index of the timestamp as of which we can start the 60sec moving window
                    dummie = 1:length(temp_t);
                    dummie = dummie(temp_t >= baseline_dt);
                    idx_baseline_dt = dummie(1); %index of the timestamp as of which we start the 10sec moving window, i.e. baseline dt
                    clear dummie temp_t
                    
                    % Index for the moving steps
                    ix1 = 1:idx_baseline_dt:length_data-idx_baseline_window; %take whole duration and then look at index you will look at (eg 1 to 60, 1 to 70 etc).
                        ix2 = idx_baseline_window:idx_baseline_dt:length_data; %index for the second bound of the window
                    ix = [ix1' ix2']; %first and last index you will check in each window (first colum onset, second colum offset) and the rows are the windows
                    if ix(end,2) < length_data
                        [ix] = [ix;[length_data-idx_baseline_window length_data]]; 
                    end
                    
                    % Calculate 8th percentile across Pre and Post datasets
                    dFF_perc8.(dFF_names{d}) = ones(1,length_data)*nan; %
                    dFF_f0.(dFF_names{d}) = ones(1,length_data)*nan;
                    
                    % calculate the entire 8th percentile (=baseline)
                    for o = 1:size(ix,1) %loop goes thru rows, equivalent of number of windows
                        idx = ix(o,1):ix(o,2); %window for this run of the for loop ; new index, saves space later. 
                        if o == 1 %
                            dFF_perc8.(dFF_names{d})(idx) = prctile(streams.ZScoredFF.(dFF_names{d})(idx),percentile); %for this first window (idx), assign 8th percentile value of dFF values during this window, 
                        else
                            dFF_perc8.(dFF_names{d})(idx(end-idx_baseline_dt)+1:end) = prctile(streams.ZScoredFF.(dFF_names{d})(idx),percentile); 
                        end
                    end
                    
                    %Substract 8th percentile, or calculated based on average to various epochs
                    %Step 1: option 1: substract moving baseline from all 
                    dFF_f0.(dFF_names{d}) = streams.ZScoredFF.(dFF_names{d}) - dFF_perc8.(dFF_names{d});

                    samples_10sec = round(10./dt_ds);
                    samples_20sec = round(20./dt_ds);
                    samples_30sec = round(30./dt_ds);
                    samples_50sec = round(50./dt_ds);
                    substract_prepost = 'moving_baseline';
                    
                    % Step 2: Option 1:
                    baseline_during1.(dFF_names{d}) = dFF_perc8.(dFF_names{d})(ix_rotarod1-samples_50sec:ix_rotarod1-samples_20sec);
                    baseline_during2.(dFF_names{d}) = dFF_perc8.(dFF_names{d})(ix_rotarod2+samples_20sec:ix_rotarod2+samples_50sec);
                    mean_baseline_during.(dFF_names{d}) = mean([baseline_during1.(dFF_names{d}) baseline_during2.(dFF_names{d})]);
                    dFF_f0.(dFF_names{d})(ix_rotarod1:ix_rotarod2+samples_20sec) = streams.ZScoredFF.(dFF_names{d})(ix_rotarod1:ix_rotarod2+samples_20sec) - mean_baseline_during.(dFF_names{d});
                    
                    %Plot
                    figure; clf
                    if show_plot == 0
                       set(gcf,'visible','off')
                    end
                    h1=plot(time_vect,streams.ZScoredFF.(dFF_names{d})-4); hold on; h2=plot(time_vect,dFF_f0.(dFF_names{d})); yline(0,'--','zero'); yline(-4,'--','zero');
                    min_x=-30; max_x=max(time_vect)+30; min_y=min(streams.ZScoredFF.(dFF_names{d}))-10; max_y=max(streams.ZScoredFF.(dFF_names{d}))+10; axis([min_x max_x min_y max_y]); 
                    xlabel('time (sec)'); ylabel('dFF'); title([dFF_names{d},' Normalized data, 8t percentile removed']);
                    h3=plot(time_vect,dFF_perc8.(dFF_names{d})-4); 

                    L=legend([h1,h2,h3],'dFF','dFF, baseline corrected','moving baseline','rotarod on/off'); L.Location = 'Best';
                    if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                        saveas(gcf,[PATH2SAVE,'figures\Zscore dFF percentilecorr ',dFF_names{d},'.tif'])
                        saveas(gcf,[PATH2SAVE,'figures\Zscore dFF percentilecorr ',dFF_names{d},'.fig'])
                    end
                    
                    %% CONVERT dFF to dFF_f0 corrected
                    streams.ZScoredFF.(dFF_names{d}) = dFF_f0.(dFF_names{d});  %% 
                    clear dFF_f0
                end
            end


            %% PEARSON R
            % Lag Analysis dFF GPe and SNr: pre/post
            Max_lag = 10; % Max time lag in seconds to test for correlation
            corr_type = 'Pearson'; % Pearson or Spearman
            save_name = 'Pearson';
              signal1 = streams.ZScoredFF.GPe(round(1/dt_ds):round((rotarodevents(1)-10)/dt_ds)); % take transpose if needed
            signal2 = streams.ZScoredFF.SNr(round(1/dt_ds):round((rotarodevents(1)-10)/dt_ds)); % take transpose if needed
            calc_plot=1;
            loop0=0;
            show_laginfo = 1;
            show_plot=1;
            save_plot=1;
            [lag2plot,Ovrl_corr] = ovrl_corr_calculation_loop...
                (signal1,signal2,time_vect,Max_lag,corr_type,calc_plot,show_plot,save_plot,PATH2SAVE,save_name,loop0,show_laginfo);
            saveas(gcf,[PATH2SAVE,'figures\dFF Pearson corr GPe SNr.fig'])
            saveas(gcf,[PATH2SAVE,'figures\dFF Pearson corr GPe SNr.tif'])
            corr_pre=max(Ovrl_corr.r,[],2);
            pearson_lag2plot_pre = Ovrl_corr.r;

            % Lag Analysis dFF GPe and SNr: during
            Max_lag = 10; % Max time lag in seconds to test for correlation
            corr_type = 'Pearson'; % Pearson or Spearman
            save_name = 'Pearson';
            signal1 = streams.ZScoredFF.GPe(round((rotarodevents(1)+60)/dt_ds):round((rotarodevents(2)-10)/dt_ds)); % take transpose if needed
            signal2 = streams.ZScoredFF.SNr(round((rotarodevents(1)+60)/dt_ds):round((rotarodevents(2)-10)/dt_ds)); % take transpose if needed
            calc_plot=1;
            loop0=0;
            show_laginfo = 1;
            show_plot=1;
            save_plot=1;
            [lag2plot,Ovrl_corr] = ovrl_corr_calculation_loop...
                (signal1,signal2,time_vect,Max_lag,corr_type,calc_plot,show_plot,save_plot,PATH2SAVE,save_name,loop0,show_laginfo);
            saveas(gcf,[PATH2SAVE,'figures\dFF Pearson corr GPe SNr.fig'])
            saveas(gcf,[PATH2SAVE,'figures\dFF Pearson corr GPe SNr.tif'])
            corr_during=max(Ovrl_corr.r,[],2);
            pearson_lag2plot_during = Ovrl_corr.r;

            % Lag Analysis dFF GPe and SNr: post
            Max_lag = 10; % Max time lag in seconds to test for correlation
            corr_type = 'Pearson'; % Pearson or Spearman
            save_name = 'Pearson';
              signal1 = streams.ZScoredFF.GPe((round(rotarodevents(2)+20)/dt_ds):(round(rotarodevents(2)+180)/dt_ds)); % take transpose if needed
            signal2 = streams.ZScoredFF.SNr((round(rotarodevents(2)+20)/dt_ds):(round(rotarodevents(2)+180)/dt_ds)); % take transpose if needed
            calc_plot=1;
            loop0=0;
            show_laginfo = 1;
            show_plot=1;
            save_plot=1;
            [lag2plot,Ovrl_corr] = ovrl_corr_calculation_loop...
                (signal1,signal2,time_vect,Max_lag,corr_type,calc_plot,show_plot,save_plot,PATH2SAVE,save_name,loop0,show_laginfo);
            saveas(gcf,[PATH2SAVE,'figures\dFF Pearson corr GPe SNr.fig'])
            saveas(gcf,[PATH2SAVE,'figures\dFF Pearson corr GPe SNr.tif'])
            corr_post=max(Ovrl_corr.r,[],2);

            corr_labels = {'pre','during','post'}
            corr = [corr_pre corr_during corr_post]
            pearson_lag2plot_post = Ovrl_corr.r;
            
            show_plot=0;
            
            %% AUC and Baseline
            % sections
            pre_start = 1; 
            pre_stop = ix_rotarod1-1;
            during_start = ix_rotarod1; 
            during_stop = ix_rotarod2-1;
            post_start = ix_rotarod2; 
            post_stop = length(time_vect);
            
            % time vector            
            time_vect_pre = time_vect(pre_start:pre_stop);
            time_vect_during = time_vect(during_start:during_stop);
            time_vect_post = time_vect(post_start:post_stop);            
            
            % GPe
            for ku=1:length(dFF_names)
                dFF.(dFF_names{ku}) = streams.ZScoredFF.(dFF_names{ku});

                % select dFF data
                dFFpre.(dFF_names{ku}) = dFF.(dFF_names{ku})(pre_start:pre_stop);
                dFFduring.(dFF_names{ku}) = dFF.(dFF_names{ku})(during_start:during_stop);
                dFFpost.(dFF_names{ku}) = dFF.(dFF_names{ku})(post_start:post_stop);


                % Calculate Baseline
                BasalLevelpre.(dFF_names{ku}) = movingbaseline(dFFpre.(dFF_names{ku}),time_vect_pre,10,5,8);  %take large window to not overestimate things - taking 10 sec instead of 30 for gcamp cos shorter duration trials
                BasalLevelduring.(dFF_names{ku}) = movingbaseline(dFFduring.(dFF_names{ku}),time_vect_during,10,5,8);
                BasalLevelpost.(dFF_names{ku}) = movingbaseline(dFFpost.(dFF_names{ku}),time_vect_post,10,5,8);

                % Find moving baseline level for each section using moving window
                MeanBasalLevelpre.(dFF_names{ku}) = mean(BasalLevelpre.(dFF_names{ku})); 
                MeanBasalLevelduring.(dFF_names{ku}) = mean(BasalLevelduring.(dFF_names{ku})); 
                MeanBasalLevelpost.(dFF_names{ku}) = mean(BasalLevelpost.(dFF_names{ku}));

                % Calculate AUC
                % Initialize min values for lower point
                % Option 1: Using minimum value. 
                Min_pre.(dFF_names{ku}) = abs(0 - min(dFFpre.(dFF_names{ku})));    
                Min_post.(dFF_names{ku}) = abs(0 - min(dFFpost.(dFF_names{ku})));
                Min_during.(dFF_names{ku}) = (Min_pre.(dFF_names{ku}) + Min_post.(dFF_names{ku}))/2;

                % Shift curve to be above min point in order to calculate AUC
                curve_pre.(dFF_names{ku}) = dFFpre.(dFF_names{ku}) + Min_pre.(dFF_names{ku});
                curve_during.(dFF_names{ku}) = dFFduring.(dFF_names{ku}) + Min_during.(dFF_names{ku});
                curve_post.(dFF_names{ku}) = dFFpost.(dFF_names{ku}) + Min_post.(dFF_names{ku});

                % Compute area under curve to lowest point
                pre_auc_uncorr.(dFF_names{ku}) = trapz(curve_pre.(dFF_names{ku})); %Trapezoidal numerical integration: Q = trapz(Y) computes the approximate integral of Y via the trapezoidal method with unit spacing
                during_auc_uncorr.(dFF_names{ku}) = trapz(curve_during.(dFF_names{ku}));
                post_auc_uncorr.(dFF_names{ku}) = trapz(curve_post.(dFF_names{ku}));

                % Normalize to epoch section duration
                pre_auc.(dFF_names{ku}) = pre_auc_uncorr.(dFF_names{ku})./(length(pre_start:pre_stop)*dt_ds);
                during_auc.(dFF_names{ku}) = during_auc_uncorr.(dFF_names{ku})./(length(during_start:during_stop)*dt_ds);
                post_auc.(dFF_names{ku}) = post_auc_uncorr.(dFF_names{ku})./(length(post_start:post_stop)*dt_ds);

                % Generate figures
                figure; clf
                if show_plot == 0
                   set(gcf,'visible','off')
                end
                s1 = subplot(3,1,1); plot(time_vect_pre,curve_pre.(dFF_names{ku}),'Color',[0.5 0.5 0.5]); hold on;
                xlabel('Time (s)'); ylabel('dFF z-score'); title(sprintf('Normalized Area under: Pre = %16.f',pre_auc.(dFF_names{ku}))); %plots the value of pre_auc: %s in scientific format; %16.f in rounded format
                % EXAMPLE: formatSpec = "The current time is: %d:%d %s"; A1 = 11; A2 = 20; A3 = 'a.m.'; str = sprintf(formatSpec,A1,A2,A3)
                area(time_vect_pre,curve_pre.(dFF_names{ku})); %area(X,Y) plots Y versus X and fills the area between 0 and Y. 
                min_x=0; max_x=max(time_vect)+30; min_y=0; max_y=max(curve_pre.(dFF_names{ku}))+10; axis([min_x max_x min_y max_y]); 
                hold off

                s2 = subplot(3,1,2); plot(time_vect_during,curve_during.(dFF_names{ku}),'Color',[0.5 0.5 0.5]);hold on
                xlabel('Time (s)'); ylabel('dFF z-score'); title(sprintf('Normalized Area under: During = %16.f',during_auc.(dFF_names{ku})))  
                area(time_vect_during,curve_during.(dFF_names{ku}));
                min_x=0; max_x=max(time_vect)+30; min_y=0; max_y=max(curve_during.(dFF_names{ku}))+10; axis([min_x max_x min_y max_y]); 
                hold off

                s3 = subplot(3,1,3); plot(time_vect_post,curve_post.(dFF_names{ku}),'Color',[0.5 0.5 0.5]); hold on
                xlabel('Time (s)'); ylabel('dFF z-score'); title(sprintf('Normalized Area under: Post = %16.f',post_auc.(dFF_names{ku})));
                area(time_vect_post,curve_post.(dFF_names{ku}));
                min_x=0; max_x=max(time_vect)+30; min_y=0; max_y=max(curve_post.(dFF_names{ku}))+10; axis([min_x max_x min_y max_y]); 
                hold off
                if done == 0 && save_plot == 1 
                    saveas(gcf,[PATH2SAVE,'figures\dFFtrial1_AUC','.tif'])
                    saveas(gcf,[PATH2SAVE,'figures\dFFtrial1_AUC','.fig'])
                end
            end
            
            
        %
        end
        
        %% SAVE INDIVDIDUAL SESSION   
        
        if BehaviorDeeplabcut == 0
            resample_Position = {};
        end
        
        if done == 0 || overwrite == 1
            Params.MICEId = mice_list(nummice).name;
            Params.SESSIONId = sessions(s).name;
            Params.Virus = Virus;
%             Params.ExperimentType = ExperimentType;
            Params.Preprocessing.lowpassfiltering = lowpassfiltering;
            Params.Preprocessing.lowpassfreq = lowpassfreq;
            Params.Preprocessing.detrending = detrending;
            Params.Preprocessing.deletemovingwindow = deletemovingwindow;
            Params.Preprocessing.highpassfiltering = highpassfiltering;
            Params.Preprocessing.highpassfreq = highpassfreq;
            Params.Path.PATH2SESSION = PATH2SESSION;
            Params.Path.PATH2SAVE = PATH2SAVE;          
            save([PATH2SAVE,'IndividualData.mat'],'Params','epocs','datatype','dFF_names',...
                'dt_ds','sampling_rate_ds','length_data','streams','time_vect',...
                'corr_pre','corr_during','corr_post','lag2plot','pearson_lag2plot_pre','pearson_lag2plot_during','pearson_lag2plot_post',...    
                'pre_auc','during_auc','post_auc','MeanBasalLevelpre','MeanBasalLevelduring','MeanBasalLevelpost'); 
        end    

    
        
        
        %% Ending all the loops

    end %for s=1:length sessions
end % for all mice: i=1:numfiles, see above

    
