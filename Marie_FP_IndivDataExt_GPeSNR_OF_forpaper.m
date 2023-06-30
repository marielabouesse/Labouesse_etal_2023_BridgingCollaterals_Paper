% Marie_FP_IndivData_extraction_GPeSNR_OF 
% DATA ANALYSIS PIPELINE FOR FIBER PHOTOMETRY EXPERIMENT: OPEN FIELD GPESNR: 
% Marie Labouesse 13 Jan 2022

% need to select the folder in which several animals data exist; in principle could have several sessions per mice but not implemented yet in detail
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
channel_names = {'GPe','SNr'}; % eg brain region
color_number = 2; % 2 if 405 and 465 recording, 3 if also 565
color_names = {'c405','c465'}; % can also be 565
dFF_number = 1; % how many dFF to analyze, can be 1 or 2, based on number of channels and colors (if more than 2, need to edit script)
dFF_names = {'GPe','SNr'}; %put them in this order: 1) channel 1, 465 and 2) channel 1, 560 OR channel 2: 465 (if other combinations, need to edit script)
datatype = {'ZScoredFF','dFF','speed'}; %data analyses of interest -- remove speed if not analysing behavior

%%%%%%%%%%% PREPROCESSING
% trimming
timetrim_start_set = 120; %time to chop off at start of recording, in seconds: adjust this manually, eg 60 seconds
timetrim_end_set = 1; %time to chop off at end of recording, eg 1 second
% low pass filtering of raw 405 and 470
lowpassfiltering = 0; % 1 if you want to low pass filter- typically yes- 0 for no
lowpassfreq = 1; % typically 1Hz
% detrending dFF
detrending = 0;% %% Write 1 for detrend dFF, normal detrend function of Matlab- typically no (Write 0) for open field- 
% substract 8th percentile as a baseline calculate using a moving window
deletemovingwindow = 0; % 1 to delete moving baseline 8th percentile, 0 if no. Typically no for open field, yes for rotarod
% high pass filtering of dFF
highpassfiltering = 0; %'yes' if you want to high pass filter, 0 for no. Typically no for open field
highpassfreq = 0.005; % typically 0.005Hz

%%%%%%%%%%% TRIAL DEFINITION
% time window for the trials and graphs
TRANGE = [-15 15]; % will create events for a -5 to +5 sec window (you can always plot less later)
% baseline correction window
BASELINE_WIN.OptoStim = [-15 -5];% baseline correction window% BASELINE_PER = [-5 -1]; 
% variables to align to
Events2Align2 = {'PulseStart_time','PulseStop_time'}; %can be {'optostim_onset','optostim_offset'} for opto, or any other TTL epoc in the "epocs" structure with same organization

%%%%%%%%%%% EXPERIMENT TYPE; NOT IN USE IN THESE OPEN FIELD ANALYSES
ExperimentType = 'power'; % Options should be 5 characters: 'power', 'freq_'. If other, then the trials are not defined and need to go in and change trial allocation.
BehaviorDeeplabcut = 1; % 0 if not, 1 if FP took camera frames itself and it will be analyzed by another program and results put into a CSV file (1 row per frame)
Opto_Powers = {'p0uW','p500uW','p2000uW'}; % POWER DEPENDENCE, dFF ; not used here
TrialType_Number = length(Opto_Powers);
stim_duration = 10; % 10 seconds or 3 seconds or 5 seconds
stim_duration2 = [];  % in case their is a ramp or another duration of interest

%%%%%%%%%%% HOW MANY MICE TO ANALYZE
LoopOrNot = 1; % 1 to loop, 0 if you only want to test one mouse and dont want to loop --> if so, edit the mouse number you want below: nummice_set
nummice_set = 1;
calccorrelations = 1;

%%%%%%%%%%% SHOW, SAVE, OVERWRITE PARAMETERS
show_plot = 0; % If 0, plots are not displayed
save_plot = 1; % If 0, plots are not saved
reanalysis = 1; % If 1, the code runs for sessions already analyzed. If 0, session excluded if analysed. This is based on the existence of the matlab space "IndividualData.mat"
overwrite = 1; % If reanalyzing data, 1 if you want to save them and overwrite results and 0 if not          
Collapse_sessions = 0; % If 0, analysis only of individual sessions. If 1, analysis of the group data (pooled),NOT IN USE HERE
Color_scale = []; % For the heatmaps. If empty, automatically adjusted to the data,NOT IN USE HERE


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
            || contains(mice_list(o).name,'data') || contains(mice_list(o).name,'figures')  || contains(mice_list(o).name,'Anymaze') ...
            || contains(mice_list(o).name,'results') || contains(mice_list(o).name,'other')  || contains(mice_list(o).name,'FILTERED')  ...
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
    AnimalID = ['ID_',mice_list(nummice).name(end-3:end)];   %last 4 digits of folders should be the animal name
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
        elseif strcmp(sessions(o).name,'data') == 1 || strcmp(sessions(o).name,'figures') == 1 || strcmp(sessions(o).name,'results') == 1 || strcmp(sessions(o).name,'Anymaze') == 1 ...
             || strcmp(sessions(o).name,'.') == 1 || strcmp(sessions(o).name,'..') == 1 || contains(sessions(o).name,'.csv') == 1 ...
            || contains(sessions(o).name,'representative figures') ||  || contains  (sessions(o).name,'FILTERED') || contains  (sessions(o).name,'other') ...
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

            %Create 'streams' structure, which will contain extracted raw data from each channel/color
            for channel=1:length(channel_names)
                for colore=1:length(color_names)
                    streams.rawdata.(channel_names{channel}).(color_names{colore}) = [];
                end
            end
            % first channel
            stream_name_type405A = {'x05A','05A','A05A'}; %possible names 
            stream_name_type465A = {'x65A','65A','A65A'}; %stick to same order on these 3 rows (names that go together on same column)
            stream_name_type560A = {'x60A','60A','A60A'};
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
                stream_name_type405B = {'x05B','05B','B05B'};
                stream_name_type465B = {'x65B','65B','B65B'};
                stream_name_type560B = {'x60A','60A','A60A'};
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
%             sampling_rate = data.streams.(stream_name_type405A{1}).fs;

            N = 10; %downsample 10 times (from initial sampling rate of 1017.3 Hz) 
            sampling_rate_ds = sampling_rate/N;
            % downsample
            for channel=1:length(channel_names)
                for colore=1:length(color_names)
                    streams.rawdata.(channel_names{channel}).(color_names{colore}) = downsample(streams.rawdata.(channel_names{channel}).(color_names{colore}),10);
                end
            end
            length_data = length(streams.rawdata.(channel_names{channel}).(color_names{colore})); %
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

            

            %% DEEPLABCUT: %MATLAB SLOW ON THIS
            % Transforming camera TTLs/frames that are beyond the end of the time vector (camera takes longer to turn off)
            if BehaviorDeeplabcut == 1 
                % last time point of recording in seconds
                time_vect_end = time_vect(length(time_vect)); % last point in time vector (in seconds)
                  frontcam_timevect_lessthan = epocs.FrontCam_start_time(epocs.FrontCam_start_time <= time_vect_end); 
                FrontCamOnLocations = zeros(1,length(frontcam_timevect_lessthan)); %time vect is already downsampled
                for i = 1:length(frontcam_timevect_lessthan)
                    if epocs.FrontCam_start_time(i) <= time_vect_end 
                        tmp = find(time_vect >= epocs.FrontCam_start_time(i));  %slow
                        FrontCamOnLocations(i) = tmp(1); 
                    end
                end

                epocs.FrontCamTTL = zeros(1,length_data);
                for i=1:length(FrontCamOnLocations)
                    epocs.FrontCamTTL(FrontCamOnLocations(i)-1) = 5; %  
                end
                clear frontcam_timevect_lessthan
                
                
                
                %% LOAD THE BEHAVIORAL INFORMATION- adjusted for DEEPLABCUT
                cameras = {'Cam1'};
                reference = {'BottomLeftCorner'};
                reference2 = {'BottomRightCorner'};
                Position =  cell(2,1);
                Body_parts =  cell(2,1);
                files_csv = dir(fullfile(PATH2SESSION,'*.csv')); 
                % Loading all the positions from both cameras
                if ~isempty(files_csv)
                    for fls = 1:length(files_csv)
                        for cam = 1:length(cameras)
                            if contains(files_csv(fls).name,cameras{cam})
                                T = readtable([files_csv(fls).folder,'\',files_csv(fls).name]);
                                for i = 2:size(T,2)
                                    tmp = table2array(T(3:end,i));
                                    Position{cam}.(T{1,i}{1}).(T{2,i}{1}) = str2double(tmp);
                                end
                                Body_parts{cam} = fieldnames(Position{cam});
                            end
                        end
                    end
                end
                % Invert y positions (camera 0,0 is on the top left and we want it at bottom left)
                for cam = 1:length(cameras)
                    tmpfld = fieldnames(Position{cam});
                    for i = 1:length(tmpfld)
                        Position{cam}.(tmpfld{i}).y = Position{cam}.(tmpfld{i}).y*-1;
                    end
                end
                % Check consistency between frame number and cam TTLs
                for cam = 1:length(cameras)
                    tmpfld = fieldnames(Position{cam});
                    tmpfld2 = fieldnames(Position{cam}.(tmpfld{1}));
                    Behav_frame_number = length(Position{cam}.(tmpfld{1}).x);
                    TTL_number = length(find(epocs.FrontCamTTL));
                    if Behav_frame_number > TTL_number
                        for i = 1:length(tmpfld)
                            for ii = 1:length(tmpfld2)
                                Position{cam}.(tmpfld{i}).(tmpfld2{ii})(TTL_number+1:end) = [];
                            end
                        end
                    end
                end      
                TTL_number_post = length(find(epocs.FrontCamTTL)) %leave visible so I can compare
                Behav_frame_number_post = length(Position{cam}.(tmpfld{1}).x) %leave visible so I can compare
 
              
                % Finding the reference points and converting positions to that reference
                for cam = 1:length(cameras)
                    tmpfld = fieldnames(Position{cam});
                    for i = 1:length(tmpfld)
                        if strcmp(tmpfld{i},reference{cam})
                            ref_val_x = median(Position{cam}.(tmpfld{i}).x(Position{cam}.(tmpfld{i}).likelihood > 0.999));
                            ref_val_y = median(Position{cam}.(tmpfld{i}).y(Position{cam}.(tmpfld{i}).likelihood > 0.999));
                        end
                    end
                    for i = 1:length(tmpfld)
                        Position{cam}.(tmpfld{i}).x = Position{cam}.(tmpfld{i}).x - ref_val_x;
                        Position{cam}.(tmpfld{i}).y = Position{cam}.(tmpfld{i}).y - ref_val_y;
                    end
                end
                % Finding the pixel to cm conversion factor
                conv_factor = cell(2,1);
                for cam = 1:length(cameras)
                    tmpfld = fieldnames(Position{cam});
                    for i = 1:length(tmpfld)
                        if strcmp(tmpfld{i},reference2{cam})
                            ref2_val_x = median(Position{cam}.(tmpfld{i}).x(Position{cam}.(tmpfld{i}).likelihood > 0.999));
                            ref2_val_y = median(Position{cam}.(tmpfld{i}).y(Position{cam}.(tmpfld{i}).likelihood > 0.999));
                        end
                    end
                    tmp_distance = sqrt(ref2_val_x^2 + ref2_val_y^2);
                    conv_factor{cam} = 42/tmp_distance; % 42cm is open field width
                    for i = 1:length(tmpfld)
                        Position{cam}.(tmpfld{i}).x = Position{cam}.(tmpfld{i}).x.*conv_factor{cam};
                        Position{cam}.(tmpfld{i}).y = Position{cam}.(tmpfld{i}).y.*conv_factor{cam};
                    end
                end
                % Calculate distance and speed for each body part
                for cam = 1:length(cameras)
                    tmpfld = fieldnames(Position{cam});
                    for i = 1:length(tmpfld)
                        Position{cam}.(tmpfld{i}).distance = ones(length(Position{cam}.(tmpfld{i}).x),1)*nan;
                        for o = 2:length(Position{cam}.(tmpfld{i}).x)
                            Position{cam}.(tmpfld{i}).distance(o) = ...
                                sqrt((Position{cam}.(tmpfld{i}).x(o) - Position{cam}.(tmpfld{i}).x(o-1))^2 + ...
                                (Position{cam}.(tmpfld{i}).y(o) - Position{cam}.(tmpfld{i}).y(o-1))^2);
                        end
                        if cam == 1
                            delta_t = diff(time_vect(FrontCamOnLocations));
                        else
                            delta_t = diff(time_vect(BackCamOnLocations));
                        end
                        delta_t = [nan;delta_t'];
                        Position{cam}.(tmpfld{i}).speed = Position{cam}.(tmpfld{i}).distance./delta_t;
                    end
                end
                
                %% Remove points with low likelihood and inpaint            
                %likelihood graph
                f= figure; f.Position = [50 70 1800 900]; sgtitle('Likelihood for each body part across time before correction')
                for cam = 1:length(cameras) %:length(cameras)
                    tmpfld = fieldnames(Position{cam});
                    for i = 1:length(tmpfld)
                        idx2extract = i;      
                        subplot(2,round(length(tmpfld)./2),i); plot(Position{cam}.(tmpfld{idx2extract}).likelihood)
                        title([tmpfld{idx2extract}])
                    end
                end
                
                %interpolation
                f = figure; f.Position = [50 70 1800 900]; sgtitle('interpolation')
                for cam = 1:length(cameras) %:length(cameras)
                    tmpfld = fieldnames(Position{cam});
                    for i = 1:length(tmpfld)
                        idx2extract = i;          
                        tmp = Position{cam}.(tmpfld{idx2extract});
                        tmp_lk = Position{cam}.(tmpfld{idx2extract}).likelihood;
                        tmp.x(tmp_lk < 0.9) = NaN;
                        tmp.x = inpaint_nans(tmp.x);
                        tmp.y(tmp_lk < 0.9) = NaN;
                        tmp.y = inpaint_nans(tmp.y);
                        tmp.likelihood(tmp_lk < 0.9) = NaN;
                        %plot (only the real body parts)
                        if idx2extract <9 %only body parts
                            subplot(2,4,i);
                            yyaxis left
                            plot(Position{cam}.(tmpfld{idx2extract}).x,'-','Color','b'); hold on
%                             plot(Position{cam}.(tmpfld{idx2extract}).y,'-','Color','b')
                            ylim([-5 45]);
                            ylabel('original data')
                            yyaxis right
                            plot(tmp.x,'-','Color','r')
%                             plot(tmp.y,'-','Color','r')
                            ylim([-5 45]);
                            ylabel('inpainted data')
                            title([tmpfld{idx2extract}])
                        end
                        %Data rewriting
%                         Position{cam}.(tmpfld{idx2extract}) = tmp;
                    end
                end
  
                %likelihood
                f= figure; f.Position = [50 70 1800 900]; sgtitle('Likelihood for each body part across time after interpolation correction (keep point if likelihood >0.9)')
                for cam = 1:length(cameras) %:length(cameras)
                    tmpfld = fieldnames(Position{cam});
                    for i = 1:length(tmpfld)
                        idx2extract = i;      
                        subplot(2,round(length(tmpfld)./2),i); plot(Position{cam}.(tmpfld{idx2extract}).likelihood)
                        title([tmpfld{idx2extract}])
                    end
                end
                
                            
                %% Plot open field trajectory of body center
                % Trajectory plot
                % new variable for animal in center or periphery zone
                cam=1;
                tmpfld = fieldnames(Position{cam});
                i=4; % body center
                Position{cam}.(tmpfld{i}).zone = NaN(length(Position{cam}.(tmpfld{i}).x),1) % center zone
                for t=1:length(Position{cam}.(tmpfld{i}).x)
                    if Position{cam}.(tmpfld{i}).x(t) >= 10 && Position{cam}.(tmpfld{i}).x(t) <= 32 && Position{cam}.(tmpfld{i}).y(t) >= 10 && Position{cam}.(tmpfld{i}).y(t) <= 32
                        Position{cam}.(tmpfld{i}).zone(t) = 0; % center zone
                    else
                        Position{cam}.(tmpfld{i}).zone(t) = 1; % periphery zone
                    end
                end
                idxcolor = Position{cam}.(tmpfld{i}).zone';

      
                % plot
                f = figure; f.Position = [100 100 500 460];
                plot(Position{cam}.(tmpfld{i}).x,Position{cam}.(tmpfld{i}).y,'Color','k','LineWidth',0.3); hold on;
                scatter(Position{cam}.(tmpfld{i}).x,Position{cam}.(tmpfld{i}).y,10,idxcolor,'Filled'); hold on;
                mycolormap = [0.3 0.7 1;1 0.3 0.1]; colormap(mycolormap)
%                 colormap(lines(2))
                min_x=0; max_x=42; min_y=0; max_y=42; axis([min_x max_x min_y max_y]); 
                xlabel('x (cm)'); ylabel('y (cm)'); title(['Body part trajectory, Animal ',AnimalID(end-3:end)]); 
                set(gca,'TickLength',[0 0])
                xticks([0 10 20 30 40]); yticks([0 10 20 30 40]);
                % open field contour
                plot([0 42],[0 0],'Color','k','LineWidth',1.5);
                plot([0 42],[42 42],'Color','k','LineWidth',1.5); 
                plot([0 0],[0 42],'Color','k','LineWidth',1.5); 
                plot([42 42],[0 42],'Color','k','LineWidth',1.5); 
                % open field center contour (10 cm)
                plot([10 32],[10 10],'Color','k','LineWidth',1.5);
                plot([10 32],[32 32],'Color','k','LineWidth',1.5); 
                plot([10 10],[10 32],'Color','k','LineWidth',1.5); 
                plot([32 32],[10 32],'Color','k','LineWidth',1.5);                 
                if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                    saveas(gcf,[PATH2SAVE,'figures\OFtrajectory.tif'])
                    saveas(gcf,[PATH2SAVE,'figures\OFtrajectory.fig'])
                end
  
                
                %% CREATE NEW STRUCTURE NewPos WITH BODY POSITION THE SAME SIZE AS FP DATA
                newPos = cell(2,1);
                for cam = 1:length(cameras)
                    if ~isempty(Position{cam})
                        tmpfld = fieldnames(Position{cam});
                        tmpfld2 = fieldnames(Position{cam}.(tmpfld{1}));
                        for i = 1:length(tmpfld)
                            for ii = 1:length(tmpfld2)
                                newPos{cam}.(tmpfld{i}).(tmpfld2{ii}) = nan(1,length_data);
                                newPos{cam}.(tmpfld{i}).(tmpfld2{ii})(epocs.FrontCamTTL(1,:) == 5) = Position{cam}.(tmpfld{i}).(tmpfld2{ii});
                            end
                        end
                    end
                end
                   
            end %end of Anymaze/Behavioral data extraction
            
            
            
            
            %% Clear TDT data structure (now everything we need is in "streams" or "epocs")
%             clear data

            %% TRIMMING           
            % Setting up trim indexes
            if BehaviorDeeplabcut == 1 %
                timetrim_start = timetrim_start_set*1; %in seconds
                timetrim_end =  timetrim_end_set*1; %in seconds
                dummie1 = 1:length(time_vect); %indexes, same number as there are samples in time_vect
                dummie2 = dummie1(time_vect > timetrim_start); %only keep the indexes starting from the new start
                idx_start = dummie2(1); %index for the new start
                dummie2 = dummie1(time_vect > time_vect(end) - timetrim_end); %only keep the indexes from the new end to the current end
                idx_end = dummie2(1); %index for the new end
                clear dummie1 dummie2 
                
            elseif BehaviorDeeplabcut == 0 %
                % trim end
                FrontCamTTL_on=find(epocs.FrontCamTTL); % 
                timetrim_end = time_vect(end)-((FrontCamTTL_on(end)+1)*dt_ds); % 
                dummie1 = 1:length(time_vect);
                dummie2 = dummie1(time_vect > time_vect(end) - timetrim_end);
                idx_end = dummie2(1); %
                timetrim_end = time_vect(idx_end); %
                timetrim_start = timetrim_start_set*1; %
                dummie2 = dummie1(time_vect > timetrim_start);
                idx_start = dummie2(1); %
                if ~isempty(Body_parts)
                    FrontCamTTL_trim = epocs.FrontCamTTL(idx_start:end);
                    FrontCamTTL_trim_on = find(FrontCamTTL_trim); %
                    idx_start = idx_start + FrontCamTTL_trim_on(1)-1; %
                    timetrim_start = time_vect(idx_start);
                end
            end
            
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
            if BehaviorDeeplabcut == 1
                for i = 1:length(tmpfld)
                    for ii = 1:length(tmpfld2)
                        newPos_trim{cam}.(tmpfld{i}).(tmpfld2{ii}) = newPos{cam}.(tmpfld{i}).(tmpfld2{ii})(:,idx_start:idx_end);
                    end
                end
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
            if BehaviorDeeplabcut == 1
                for cam = 1:length(cameras)
                    if ~isempty(newPos_trim{cam})
                        % Ytail = Ytail_trim;
                        Position{cam} = newPos_trim{cam};
                    end
                end
                clear newPos_trim
            end
            length_data = length(streams.rawdata.(channel_names{channel}).(color_names{colore}));
            clear epocs_trim time_vect_trim 
            streams = rmfield(streams,'rawdata_trim'); %clear variable within structure.

            
            % Need to do FrontCamTTL again after trimming
                time_vect_end = time_vect(length(time_vect)); 
                frontcam_timevect_lessthan = epocs.FrontCam_start_time(epocs.FrontCam_start_time <= time_vect_end); 
                FrontCamOnLocations = zeros(1,length(frontcam_timevect_lessthan)); 
                for i = 1:length(frontcam_timevect_lessthan)
                    if epocs.FrontCam_start_time(i) <= time_vect_end 
                        tmp = find(time_vect >= epocs.FrontCam_start_time(i));  %slow
                        FrontCamOnLocations(i) = tmp(1); 
                    end
                end

                epocs.FrontCamTTL = zeros(1,length(time_vect));
                for i=1:length(FrontCamOnLocations)
                    epocs.FrontCamTTL(FrontCamOnLocations(i)) = 5;    
                end
                clear frontcam_timevect_lessthan

                
                
            %% Upsample the body position (Anymaze)
            if BehaviorDeeplabcut == 1
                resample_Position = cell(2,1);
                for cam = 1:length(cameras)
                    if ~isempty(Position{cam})
                        tmpfld = fieldnames(Position{cam});
                        tmpfld2 = fieldnames(Position{cam}.(tmpfld{1}));
                        for i = 1:length(tmpfld)
                            for ii = 1:length(tmpfld2)
                                resample_Position{cam}.(tmpfld{i}).(tmpfld2{ii}) = inpaint_nans(Position{cam}.(tmpfld{i}).(tmpfld2{ii}));
                            end
                        end
                    end
                end
            end

            
            %% Extract specific body part position: SPEED
            if BehaviorDeeplabcut == 1
                if ~isempty(Body_parts)
                    cam = 1;
                    BodySpeed = resample_Position{cam}.BodyCenter.speed;  %

                    figure; clf
                    if show_plot == 0
                           set(gcf,'visible','off')
                    end
                    plot(time_vect,Position{cam}.BodyCenter.speed,'.b') % raw data from Anymaze 
                    hold on
                    plot(time_vect,BodySpeed,'r') % upsampled data
                    yline(0.03)
                    yline(0.06)
                    min_x=-30; max_x=max(time_vect)+30; min_y=min(BodySpeed)-1; max_y=max(BodySpeed)+1; axis([min_x max_x min_y max_y]); 
                    xlabel('time(sec)'); ylabel('speed (cm/s)'); title('Speed data'); 
                    L=legend('Speed raw','Speed upsampled'); L.Location = 'Best';
                %     BehavData = resample_Position; %put in BehavData the upsampled data

                    if  done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                        saveas(gcf,[PATH2SAVE,'figures\Speed upsampling, ',AnimalID,'.tif'])
                        saveas(gcf,[PATH2SAVE,'figures\Speed upsampling, ',AnimalID,'.fig'])
                    end
                end
            end
            % clear resample_Position field_names

            %% LOW PASS FILTER OF FP DATA TO 1Hz
            if lowpassfiltering == 1
                %filter option 2: butter
                ftype = 'low';
                n_trials = 2; % 2nd order filter, like Arturo and Thomas Akan
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
            % Initialize
            F0_405 = []; F0_465 = []; F0_560 = []; F0_405_2 = []; F0_465_2 = [];
            dFF_405 = []; dFF_465 = []; dFF_560 = []; dFF_405_2 = []; dFF_465_2 = [];
            % Calculations in each channel and color
            % First channel, 465 color, first dFF
            channel=1;
                % % dFF synapse; % Note: Synapse uses an exponential smooth to estimate the F0 (I dont)
            F0_405 = movingaverage(streams.rawdata.(channel_names{channel}).c405,time_vect,120,10); % sliding average in a moving window
            F0_465 = movingaverage(streams.rawdata.(channel_names{channel}).c465,time_vect,120,10); % sliding average in a moving window
            % channel dFF calculation
            dFF_405 = (streams.rawdata.(channel_names{channel}).c405-F0_405)./F0_405;
            dFF_465 = (streams.rawdata.(channel_names{channel}).c465-F0_465)./F0_465;
            % final dFF calculation
            streams.dFFsy.(dFF_names{1}) = 100*(dFF_465-dFF_405);
                % % dFF polyfit
            % Step 1: Polyfit on entire data 
            index = 1:length(streams.rawdata.(channel_names{channel}).c465);
            calc_coeff_fit_465 = polyfit(streams.rawdata.(channel_names{channel}).c405,streams.rawdata.(channel_names{channel}).c465,1); %fitting control Dv2 against signal Dv1 using least squares fit of degree 1
            % Step 2: Apply polyfit on 465 
            data405_fitted_465 = calc_coeff_fit_465(1).*streams.rawdata.(channel_names{channel}).c405 + calc_coeff_fit_465(2);
            % Step 3: Normalize
            streams.dFF.(dFF_names{1}) = 100*((streams.rawdata.(channel_names{channel}).c465 - data405_fitted_465)./data405_fitted_465); % of deltaF/F              

            % second dFF, either 2nd color in channel 1 or 2nd channel (1 color)- if ever I had 2 color, 2 channels, would need to edit
            % this script
            if color_number == 3 && channel_number == 1
                channel = 1;
                    % % dFF synapse; % Note: Synapse uses an exponential smooth to estimate the F0 (I dont)
                F0_405 = movingaverage(streams.rawdata.(channel_names{channel}).c405,time_vect,120,10); % sliding average in a moving window
                F0_560 = movingaverage(streams.rawdata.(channel_names{channel}).c560,time_vect,120,10); % sliding average in a moving window
                % channel dFF calculation
                dFF_405 = (streams.rawdata.(channel_names{channel}).c405-F0_405)./F0_405;
                dFF_560 = (streams.rawdata.(channel_names{channel}).c560-F0_560)./F0_560;
                % final dFF calculation
                streams.dFFsy.(dFF_names{2}) = 100*(dFF_560-dFF_405);
                    % % dFF polyfit
                % Step 1: Polyfit on entire data
                index = 1:length(streams.rawdata.(channel_names{channel}).c560);
                calc_coeff_fit_560 = polyfit(streams.rawdata.(channel_names{channel}).c405,streams.rawdata.(channel_names{channel}).c560,1); %fitting control Dv2 against signal Dv1 using least squares fit of degree 1
                % Step 2: Apply polyfit on 560 
                data405_fitted_560 = calc_coeff_fit_560(1).*streams.rawdata.(channel_names{channel}).c405 + calc_coeff_fit_560(2);
                % Step 3: Normalize
                streams.dFF.(dFF_names{2}) = 100*((streams.rawdata.(channel_names{channel}).c560 - data405_fitted_560)./data405_fitted_560); % of deltaF/F 
            elseif color_number == 2 && channel_number == 2
                channel = 2;
                    % % dFF synapse; % Note: Synapse uses an exponential smooth to estimate the F0 (I dont)
                F0_405_2 = movingaverage(streams.rawdata.(channel_names{channel}).c405,time_vect,120,10); % sliding average in a moving window
                F0_465_2 = movingaverage(streams.rawdata.(channel_names{channel}).c465,time_vect,120,10); % sliding average in a moving window
                % channel dFF calculation
                dFF_465_2 = (streams.rawdata.(channel_names{channel}).c465-F0_465_2)./F0_465_2;
                dFF_405_2 = (streams.rawdata.(channel_names{channel}).c405-F0_405_2)./F0_405_2;
                % final dFF calculation
                streams.dFFsy.(dFF_names{2}) = 100*(dFF_465_2-dFF_405_2);
                    % % dFF polyfit
                % Step 1: Polyfit on entire data 
                index = 1:length(streams.rawdata.(channel_names{channel}).c465);
                calc_coeff_fit_465_2 = polyfit(streams.rawdata.(channel_names{channel}).c405,streams.rawdata.(channel_names{channel}).c465,1); %fitting control Dv2 against signal Dv1 using least squares fit of degree 1
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
                plot(time_vect,streams.dFFsy.(dFF_names{d})); 
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

            %% ADJUSTING dFF BASELINE TO 0 IN A MOVING WINDOW      
            if deletemovingwindow == 1;
                for d=1:length(dFF_names)
                    % Define parameters for moving window 
                    baseline_window = 10; %adjust 
                    baseline_dt = 5; %Moving steps for the threshold 
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
                    ix1 = 1:idx_baseline_dt:length_data-idx_baseline_window; 
                    ix2 = idx_baseline_window:idx_baseline_dt:length_data; %index for the second bound of the window
                    ix = [ix1' ix2']; %first and last index you will check in each window (first colum onset, second colum offset) and the rows are the windows
                    if ix(end,2) < length_data
                        [ix] = [ix;[length_data-idx_baseline_window length_data]]; %
                    end
                    % Calculate 8th percentile across Pre and Post datasets
                    dFF_perc8.(dFF_names{d}) = ones(1,length_data)*nan; 
                    dFF_f0.(dFF_names{d}) = ones(1,length_data)*nan;
                    % calculate the entire 8th percentile (=baseline)
                    for o = 1:size(ix,1) %
                        idx = ix(o,1):ix(o,2); %
                        if o == 1
                            dFF_perc8.(dFF_names{d})(idx) = prctile(streams.dFF.(dFF_names{d})(idx),percentile); %for this first window (idx), assign 8th percentile value of dFF values during this window, 
                        else
                            dFF_perc8.(dFF_names{d})(idx(end-idx_baseline_dt)+1:end) = prctile(streams.dFF.(dFF_names{d})(idx),percentile); 
                        end
                    end
                    %Substract 8th percentile, or calculated based on average to various epochs
                    %Step 1: option 1: substract moving baseline from all 
                    dFF_f0.(dFF_names{d}) = streams.dFF.(dFF_names{d}) - dFF_perc8.(dFF_names{d});

                    %Plot
                    figure; clf
                    if show_plot == 0
                       set(gcf,'visible','off')
                    end
                    h1=plot(time_vect,streams.dFF.(dFF_names{d})-4); hold on; h2=plot(time_vect,dFF_f0.(dFF_names{d})); yline(0,'--','zero'); yline(-4,'--','zero');
                    min_x=-30; max_x=max(time_vect)+30; min_y=min(streams.dFF.(dFF_names{d}))-10; max_y=max(streams.dFF.(dFF_names{d}))+10; axis([min_x max_x min_y max_y]); 
                    xlabel('time (sec)'); ylabel('dFF'); title([dFF_names{d},' Normalized data, 8t percentile removed']);
                    h3=plot(time_vect,dFF_perc8.(dFF_names{d})-4); 
                    for o = 1:size(epocs.PulseStart_time)
                        plot([epocs.PulseStart_time(o) epocs.PulseStart_time(o)],[-5 5],'m','linewidth',1)
                        plot([epocs.PulseStop_time(o) epocs.PulseStop_time(o)],[-5 5],'c','linewidth',1)
                    end
                    L=legend([h1,h2,h3],'dFF','dFF, baseline corrected','moving baseline','rotarod on/off'); L.Location = 'Best';
                    if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                        saveas(gcf,[PATH2SAVE,'figures\dFF percentilecorr ',dFF_names{d},'.tif'])
                        saveas(gcf,[PATH2SAVE,'figures\dFF percentilecorr ',dFF_names{d},'.fig'])
                    end
                    %% CONVERT dFF to dFF_f0 corrected
                    streams.dFF.(dFF_names{d}) = dFF_f0.(dFF_names{d});
                    clear dFF_f0
                end
            end

            
            %% HIGH-PASS FILTERING
            if highpassfiltering == 1; %'yes' if you want to high pass filter- typically off for open field
                % Filter out the slow fluctuation
                ftype = 'high';
                n_trials = 2; % 2nd order filter, like Arturo and Thomas Akan
                Wn = highpassfreq/((sampling_rate_ds)/2); %highpassfreq:defined above
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
            %% Calculations
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

            %% ONCE DONE WITH MODIFYNG dFF in streams, add SPEED (can also modify later)          
            if BehaviorDeeplabcut == 1
                for d=1:length(dFF_names)   
                    streams.speed.(dFF_names{d}) = BodySpeed;
                end
            end
            
                       
            
            %% Correlations between dFF and Behavioral data                     
          
            if calccorrelations == 1
                for d=1:length(dFF_names)                
%                     BodySpeed_cms = BodySpeed*100;
                    % Plot raw data version dFFs and speed together
                    f=figure; clf; subplot(1,2,1); f.Position = [100 100 1200 800];
                    if show_plot == 0
                       set(gcf,'visible','off')
                    end
                    yyaxis left
                    h1=plot(time_vect,streams.dFF.(dFF_names{d})); hold on;
                    ylabel('dFF'); 
                    min_x=-30; max_x=max(time_vect)+30; min_y=min(streams.dFF.(dFF_names{d}))-2; max_y=max(streams.dFF.(dFF_names{d}))+2; axis([min_x max_x min_y max_y]); 
                    for o = 1:size(epocs.PulseStart_time)
                        plot([epocs.PulseStart_time(o) epocs.PulseStart_time(o)],[-1 1],'-','color','k','linewidth',1)
                        plot([epocs.PulseStop_time(o) epocs.PulseStop_time(o)],[-1 1],'-','color','k','linewidth',1)
                    end
                    yyaxis right
                    h2=plot(time_vect,BodySpeed);
                    ylabel('Speed (cm/s)'); 
                    yline(0,'--','zero'); 
                    min_x=-30; max_x=max(time_vect)+30; min_y=min(BodySpeed)-5; max_y=max(BodySpeed)+5; axis([min_x max_x min_y max_y]); 
                    xlabel('time (sec)'); title([dFF_names{d},' dFF and Mouse Speed']);
                    L=legend([h1,h2],'dFF','Speed'); L.Location = 'Best';
                    if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                        saveas(gcf,[PATH2SAVE,'figures\dFF and speed ',dFF_names{d},'.tif'])
                        saveas(gcf,[PATH2SAVE,'figures\dFF and speed ',dFF_names{d},'.fig'])
                    end
                    % Plot smoothed versions dFFs and speed together
                    smooth_window2 = 1./dt_ds; % first element here is in seconds
                    dFF_smooth.(dFF_names{d}) = streams.dFF.(dFF_names{d});
                    BodySpeed_smooth = smooth(BodySpeed,smooth_window2);
                    %plot
                    subplot(1,2,2);
                    if show_plot == 0
                       set(gcf,'visible','off')
                    end
                    yyaxis left
                    h1=plot(time_vect,dFF_smooth.(dFF_names{d})); hold on;
                    ylabel('dFF'); 
                    min_x=-30; max_x=max(time_vect)+30; min_y=min(dFF_smooth.(dFF_names{d}))-2; max_y=max(dFF_smooth.(dFF_names{d}))+2; axis([min_x max_x min_y max_y]); 

                    yyaxis right
                    h2=plot(time_vect,BodySpeed_smooth);
                    ylabel('Speed (cm/s)'); 
                    yline(0,'--','zero'); 
                    min_x=-30; max_x=max(time_vect)+30; min_y=min(BodySpeed_smooth)-5; max_y=max(BodySpeed_smooth)+5; axis([min_x max_x min_y max_y]); 
                    xlabel('time (sec)'); title([dFF_names{d},'dFF and Mouse Speed, smoothed']);
                    L=legend([h1,h2],'dFF','Speed smoothed'); L.Location = 'Best';
                    smooth_window_fortitle = sprintf('%.0f',smooth_window2*dt_ds)
                    title(['same with speed smoothed (',smooth_window_fortitle,' s)'])
                    if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                        saveas(gcf,[PATH2SAVE,'figures\dFF and speed smoothed ',dFF_names{d},'.fig'])
                        saveas(gcf,[PATH2SAVE,'figures\dFF and speed smoothed ',dFF_names{d},'.tif'])
                    end
                    
                    
                    %% Lag Analysis dFF and speed at different speed smoothing 
                    smoothwindows = [0.1 0.2 0.5 1 2 5 10 20 50]; % in seconds
                    allcorr = NaN(1,length(smoothwindows));
                    for z= 1:length(smoothwindows)
                        smooth_window = smoothwindows(z)./dt_ds; % first element here is in seconds
                        smooth_window_sec = sprintf('%.0f',smooth_window*dt_ds)
                        dFF_smooth.(dFF_names{d}) = streams.dFF.(dFF_names{d});
                        BodySpeed_smooth = smooth(BodySpeed,smooth_window);
                        
                        Max_lag = 10; % Max time lag in seconds to test for correlation
                        corr_type = 'Pearson'; % Pearson or Spearman
                        save_name = 'Pearson';
                        signal1 = streams.dFF.(dFF_names{d}); % take transpose if needed
                        signal2 = BodySpeed_smooth'; % take transpose if needed
                        if show_plot == 0
                           set(gcf,'visible','off')
                        end
                        calc_plot=1;
                        loop0=0;
                        show_laginfo = 1;
                        [lag2plot,Ovrl_corr] = ovrl_corr_calculation_loop...
                            (signal1,signal2,time_vect,Max_lag,corr_type,calc_plot,show_plot,save_plot,PATH2SAVE,save_name,loop0,show_laginfo);
                        if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                            saveas(gcf,[PATH2SAVE,'figures\dFF speed Pearson corr ',smooth_window_sec,'s smooth ',dFF_names{d},'.tif'])
                            saveas(gcf,[PATH2SAVE,'figures\dFF speed Pearson corr ',smooth_window_sec,'s smooth',dFF_names{d},'.tif'])
                        end
                        allcorr(z)=max(Ovrl_corr.r,[],2);
                    end
                    % average plot
                    f=figure;
                    plot(smoothwindows,allcorr,'Color','k','Linewidth',1,'Linestyle','-'); hold on; scatter(smoothwindows,allcorr,'Filled');
                    min_x=-1; max_x=60; min_y=0.2; max_y=0.6; axis([min_x max_x min_y max_y]); 
                    xlabel('Speed smoothing bin size (sec)'); ylabel('Correlation dFF mouse speed (r)'); 
                    title([dFF_names{d},' dFF and speed correlation at different speed smooth bins']);
                    txt = ['max corr = ',num2str(max(allcorr,[],2)),', at ',num2str(smoothwindows(find(allcorr==max(allcorr,[],2)))),'sec bins'];
%                     text(80,80,txt,'HorizontalAlignment','right')
                    L=legend([txt]); L.Location = 'Best';
                    if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                            saveas(gcf,[PATH2SAVE,'figures\dFF speed Pearson corr overall',dFF_names{d},'.fig'])
                            saveas(gcf,[PATH2SAVE,'figures\dFF speed Pearson corr overall',dFF_names{d},'.tif'])
                    end                    

              
                    
                end
            
         
            
                 %% Lag Analysis dFF GPe and SNr 
                Max_lag = 10; % Max time lag in seconds to test for correlation
                corr_type = 'Pearson'; % Pearson or Spearman
                save_name = 'Pearson';
                signal1 = streams.ZScoredFF.GPe; % take transpose if needed
                signal2 = streams.ZScoredFF.SNr; % take transpose if needed
                calc_plot=1;
                loop0=0;
                show_laginfo = 1;
                if show_plot == 0
                   set(gcf,'visible','off')
                end
                [lag2plot,Ovrl_corr] = ovrl_corr_calculation_loop...
                    (signal1,signal2,time_vect,Max_lag,corr_type,calc_plot,show_plot,save_plot,PATH2SAVE,save_name,loop0,show_laginfo);
                if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                    saveas(gcf,[PATH2SAVE,'figures\dFF Pearson corr GPe SNr.fig'])
                    saveas(gcf,[PATH2SAVE,'figures\dFF Pearson corr GPe SNr.tif'])
                end
                allcorr(z)=max(Ovrl_corr.r,[],2);

            end
            
            
            
            
            %% TRIGGER MOVEMENT ONSET AND OFFSET
            
            %% FIND MOVEMENT BOUTS OF THE BODY CENTER
            z=2; % 2 sec
            smooth_window = z./dt_ds;
            BodySpeed_smooth = smooth(BodySpeed,smooth_window);
            speed4analysis = BodySpeed_smooth'; % 
            % Define events based on threshold crossing, ie which are above the first threshold
            thshl1 = 3; % In cm/s, find bouts above the threshold
            thshl2 = 5; % The event must surpass this value at least once to keep it
            tmp = speed4analysis >= thshl1; % 1 if body speed above 2cm, 0 otherwise
            A = find(tmp == 1); % index of when body speed above 2
            B = diff(A); % 
            C = find(B ~= 1); % find points where B is not equal to 1, ie the changing points 
            onset_mov = A(C+1); % fixnd points in A with the index C (changing points) + 1 cos that's the onset
            onset_mov = [A(1), onset_mov]; % also take the first point in A which is the first onset
            offset_mov = A(C); % find points in A with the index C which are all the offsets
            offset_mov = [offset_mov,A(length(A))]; % also take the last point in A which is the last offset
            figure; scatter([1:length(offset_mov)],(offset_mov-onset_mov)*dt_ds); ylabel('movement duration(s)'); xlabel('movement unit')
%             
%             % delete the events where speed never reaches the second threshold
            for o = length(onset_mov):-1:1
                if sum(speed4analysis(onset_mov(o):offset_mov(o)) > thshl2) == 0 % i.e. you take the sum to see if you reach at least once a 1 value above the threshold
                    onset_mov(o) = [];
                    offset_mov(o) = [];
                else 
                    tmp2 = find(speed4analysis(onset_mov(o):offset_mov(o)) > thshl2);
                    tmp3 = length(speed4analysis(onset_mov(o):offset_mov(o)));
                    offset_mov(o) = offset_mov(o) - (tmp3 - tmp2(end));     
                end
            end
            figure; scatter([1:length(offset_mov)],(offset_mov-onset_mov)*dt_ds); ylabel('movement duration(s)'); xlabel('movement unit')
            
            % identify the onset of movements
            show_plot = 1;
            onset_mov_real = ones(length(onset_mov),1)*nan;
            lookbacktime = 350; % 300 frames is 3 sec
            figure; hold on;
            for o = 1:length(onset_mov)
                if onset_mov(o)-lookbacktime>0
                    temp = speed4analysis(onset_mov(o)-lookbacktime:onset_mov(o)+lookbacktime); % look back 3 sec; at 100fps, that is 300
                    temp = -temp; % turn it around so that findpeaks allows to find valleys
                   time_vect_temp = time_vect(onset_mov(o)-lookbacktime:onset_mov(o)+lookbacktime);
                   [pks,idx] = findpeaks(temp,time_vect_temp,'MinPeakProminence',0.3);
                    if length(pks)>1
                        pks2=-min(-pks);
                        idx2=find(pks == pks2); 
                        idx = idx(idx2);
                        pks=pks2;
                    elseif length(pks) == 0
                       pks = min(temp);
                       idx = time_vect_temp(find(temp==pks));
                    end
                    index_timetemp = find(time_vect_temp == idx); %index in time_vect_temp where we got our point
                    onset_mov_real(o) = onset_mov(o)-lookbacktime+index_timetemp;
                    if idx > time_vect(onset_mov(o))
                        idx = time_vect(onset_mov(o));
                        onset_mov_real(o) = onset_mov(o);
                    end
                    if speed4analysis(onset_mov_real(o))>speed4analysis(onset_mov(o))
                        idx = time_vect(onset_mov(o));
                        onset_mov_real(o) = onset_mov(o);
                    end                    
                    plot(time_vect(onset_mov(o)-lookbacktime:onset_mov(o)+lookbacktime),-temp); hold on; 
    %                 plot(c,-pks,'*c','LineWidth',4);
                    plot(time_vect(onset_mov(o)),speed4analysis(onset_mov(o)),'*g','LineWidth',3)
                    plot(time_vect(onset_mov_real(o)),speed4analysis(onset_mov_real(o)),'*c','LineWidth',3)
                    clear a b c
                end
            end
            % test
            onset_mov = onset_mov_real;
           
            
             % identify the offset of movements
            show_plot = 1;
            offset_mov_real = ones(length(onset_mov),1)*nan;

            figure; hold on;
            for o = 1:length(offset_mov)
                lookbacktime = 300; % 300 frames is 3 sec
                if lookbacktime > offset_mov(o)-onset_mov(o)
                    lookbacktime = offset_mov(o)-onset_mov(o)
                end
                lookaheadtime = 100;
                if offset_mov(o)+lookaheadtime > length(time_vect)
                    lookaheadtime = 0;
                end
                if onset_mov(o)-lookbacktime>0
                    temp = speed4analysis(offset_mov(o)-lookbacktime:offset_mov(o)+lookaheadtime); % look back 3 sec; at 100fps, that is 300
                    time_vect_temp = time_vect(offset_mov(o)-lookbacktime:offset_mov(o)+lookaheadtime);
                    temp2 = max(temp);
                    idx_temp = find(temp == temp2);
                       offset_mov_real(o) = offset_mov(o)-lookbacktime+idx_temp;
                    if onset_mov(o)>offset_mov_real(o)
                        lookbacktime = offset_mov(o)-onset_mov(o); % 300 frames is 3 sec
                        lookaheadtime = 100;
                        if offset_mov(o)+lookaheadtime > length(time_vect)
                            lookaheadtime = 0;
                        end
                        temp = speed4analysis(offset_mov(o)-lookbacktime:offset_mov(o)+lookaheadtime); % look back 3 sec; at 100fps, that is 300
                        time_vect_temp = time_vect(offset_mov(o)-lookbacktime:offset_mov(o)+lookaheadtime);
                        temp2 = max(temp);
                        idx_temp = find(temp == temp2);
                    end
                       if offset_mov(o)+lookaheadtime < length(time_vect)
                        plot(time_vect(offset_mov(o)-lookbacktime:offset_mov(o)+lookaheadtime),temp); hold on; 
        %                 plot(c,-pks,'*c','LineWidth',4);
                        plot(time_vect(offset_mov(o)),speed4analysis(offset_mov(o)),'*g','LineWidth',3)
                        plot(time_vect(offset_mov_real(o)),speed4analysis(offset_mov_real(o)),'*c','LineWidth',3)
                        clear a b c
                    end
                end
            end
            
            
            % corret onset and offset
            for o = length(onset_mov):-1:1
                if isnan(onset_mov(o))
                    onset_mov(o) = [];
                    offset_mov_real(o) = [];
                elseif isnan(offset_mov_real(o))
                    onset_mov(o) = [];
                    offset_mov_real(o) = [];
                end
            end
             % test
            figure; 
            plot(time_vect,speed4analysis); hold on;
            plot(time_vect(offset_mov_real),speed4analysis(offset_mov_real),'*b','LineWidth',3)
            plot(time_vect(offset_mov),speed4analysis(offset_mov),'*g','LineWidth',3)
            plot(time_vect(onset_mov),speed4analysis(onset_mov),'*r','LineWidth',3)
            yline(thshl1,':m','LineWidth',2);
            yline(thshl2,'-.m','LineWidth',2);
            % convert
            offset_mov = offset_mov_real;
            for o = 1:length(onset_mov)
                if onset_mov(o)>offset_mov(o)
                    display('error in the calculation of onset and offset, offset has a smaller index than onset')                 
                end
            end
            
            
 

            % define vectors that give, for each indidivual movement, the maximal speed, the time vector and the auc
            peak_speed_mov = ones(length(onset_mov),1)*nan;
            dur_speed_mov = ones(length(onset_mov),1)*nan;
            auc_speed_mov = ones(length(onset_mov),1)*nan;
            for o = 1:length(onset_mov)
                peak_speed_mov(o) = max(speed4analysis(onset_mov(o):offset_mov(o)));
                dur_speed_mov(o) = time_vect(length(onset_mov(o):offset_mov(o)));
                auc_speed_mov(o) = trapz(speed4analysis(onset_mov(o):offset_mov(o)));
            end

            body_mov.onset = onset_mov;
            body_mov.offset = offset_mov;
            body_mov.peak = peak_speed_mov;
            body_mov.dur = dur_speed_mov;
            body_mov.auc = auc_speed_mov;

            %% plot the speed over time with identification of the movement onset and offset
            figure
            if show_plot == 0
                set(gcf,'visible','off')
            end
            subplot(2,1,1)            
            plot(time_vect,BodySpeed,'k','LineStyle','-','LineWidth',1.5); hold on
            plot(time_vect,BodySpeed_smooth,'b','LineStyle','-','LineWidth',1.5); hold on
            yline(thshl1,':m','LineWidth',2);
            yline(thshl2,'-.m','LineWidth',2);
            xlim([1 time_vect(end-1)])
            xlabel('Time (s)')
            ylabel('Mouse speed (cm/s)')
            hold on
            plot(time_vect(onset_mov),speed4analysis(onset_mov),'*g','LineWidth',3)
            plot(time_vect(offset_mov),speed4analysis(offset_mov),'*r','LineWidth',3)
            title('Selected movements')
            legend('Mouse speed','Mouse speed smoothed','Threshold 1','Threshold 2','Movement onset',...
                'Movement offset','dFF GPe','dFF SNr','Location','Best')
            axis([100 200 -4 20]);            
            
            subplot(2,1,2) 
            yyaxis left
            plot(time_vect,BodySpeed_smooth,'b','LineStyle','-','LineWidth',1.5); hold on
            ylabel('Mouse speed (cm/s)')
            axis([0 800 -4 20]);            

            yyaxis right
            ylabel('Zscore dFF')
            plot(time_vect,streams.ZScoredFF.GPe,'LineStyle','-','Color',[0.4940 0.1840 0.5560],'LineWidth',1.5); hold on
            plot(time_vect,streams.ZScoredFF.SNr,'LineStyle','-','Color',[0.4660 0.6740 0.1880],'LineWidth',1.5)    
            axis([100 200 -3 3]);            
            if save_plot
                saveas(gcf,[PATH2SAVE,'figures\Selected movements'],'tif');
                saveas(gcf,[PATH2SAVE,'figures\Selected movements'],'fig');
            end

            
        
            %% Align the dFF to head movements: movement onset or offset period
            TRANGE = [-5 5];
            BASELINE_PER = [-0.2 0];
            dummie = 1:length(time_vect);
            dummie = dummie(time_vect >= abs(TRANGE(1)));
            idx_Init = dummie(1); % index for first column of trial, will be 1 if trial starts at 0, will be negative if trials start before 0
            dummie = 1:length(time_vect);
            dummie = dummie(time_vect >= abs(TRANGE(2)));
            idx_End = dummie(1); % last index for trial, eg 500 for a 5 sec trial at 100 fps

            n_trials = idx_Init + idx_End; % Length of each trial
            t_trials = time_vect(1:n_trials) - time_vect(idx_Init);

            Events = {'mov_onset','mov_offset'};
            for i = 1:length(Events)
                Stim_data.(Events{i}).idx = ones(size(body_mov.onset,2),n_trials)*nan;
                % Stim_data.tail_jump.corr = ones(size(valley.loc,1),n)*nan;
                Stim_data.(Events{i}).dFF.(dFF_names{1}).raw = ones(size(body_mov.onset,2),n_trials)*nan;
                Stim_data.(Events{i}).dFF.(dFF_names{1}).baseline_corrected = ones(size(body_mov.onset,2),n_trials)*nan;
                Stim_data.(Events{i}).dFF.(dFF_names{2}).raw = ones(size(body_mov.onset,2),n_trials)*nan;
                Stim_data.(Events{i}).dFF.(dFF_names{2}).baseline_corrected = ones(size(body_mov.onset,2),n_trials)*nan;
                Stim_data.(Events{i}).body_speed.raw = ones(size(body_mov.onset,2),n_trials)*nan;
                Stim_data.(Events{i}).body_speed.baseline_corrected = ones(size(body_mov.onset,2),n_trials)*nan;
                Stim_data.(Events{i}).firstderiv.(dFF_names{1}).raw = ones(size(body_mov.onset,2),n_trials-1)*nan; % 1 column less
                Stim_data.(Events{i}).firstderiv.(dFF_names{1}).baseline_corrected = ones(size(body_mov.onset,2),n_trials-1)*nan; % 1 column less
                Stim_data.(Events{i}).firstderiv.(dFF_names{2}).raw = ones(size(body_mov.onset,2),n_trials-1)*nan; % 1 column less
                Stim_data.(Events{i}).firstderiv.(dFF_names{2}).baseline_corrected = ones(size(body_mov.onset,2),n_trials-1)*nan; % 1 column less
                
            end
            for i = 1:length(Events) % onset and offset
                for o = 1:size(body_mov.onset,1)
                    if i == 1 % movement onset
                        ix = body_mov.onset(o); % follow the indices of trials in body_mov.onset (index in the time_vect)
                    else
                        ix = body_mov.offset(o); % follow the indices of trials in body_mov.offset (index in the time_vect)
                    end
                    if ix+idx_End <= length(time_vect) && ix-(idx_Init-1) >= 1 % for trials within the boundaries of time_vect
                        tmp = ix - (idx_Init-1):ix + idx_End; %indices for the size of trials we want, indices correspond to indices in time_vect
                        Stim_data.(Events{i}).idx(o,:) = tmp;
                        for ii = 1:channel_number % GPe and SNr
                            Stim_data.(Events{i}).dFF.(dFF_names{ii}).raw(o,:) = ...
                                streams.ZScoredFF.(dFF_names{ii})(tmp);
                            % first derivative:
                            dy = diff(Stim_data.(Events{i}).dFF.(dFF_names{ii}).raw(o,:)); % differences between contiguous elements: substract point n to the point in n+1 --> 1 point smaller than original data
                            dx = diff(time_vect(tmp)); % differences between contiguous elements, time vector. doesnt matter if starts at 0. will have 1 element less                           
                            Stim_data.(Events{i}).firstderiv.(dFF_names{ii}).raw(o,1:length(tmp)-1) = dy./dx; % first derivative  - last point is a NaN
                        end
                        Stim_data.(Events{i}).body_speed.raw(o,:) = speed4analysis(tmp);
                    elseif ix+idx_End > length(time_vect) % last trials
                        tmp = ix - (idx_Init-1):length(time_vect)-1;
                        Stim_data.(Events{i}).idx(o,1:length(tmp)) = tmp;
                        for ii = 1:length(channel_number)
                            Stim_data.(Events{i}).dFF.(dFF_names{ii}).raw(o,1:length(tmp)) = ...
                                streams.ZScoredFF.(dFF_names{ii})(tmp);
                            % first derivative:
                            dy = diff(Stim_data.(Events{i}).dFF.(dFF_names{ii}).raw(o,1:length(tmp))); % differences between contiguous elements: substract point n to the point in n+1 --> 1 point smaller than original data
                            dx = diff(time_vect(tmp)); % differences between contiguous elements, time vector. doesnt matter if starts at 0. will have 1 element less                           
                            Stim_data.(Events{i}).firstderiv.(dFF_names{ii}).raw(o,1:length(tmp)-1) = dy./dx; % first derivative  - last point is a NaN
                        end
                        Stim_data.(Events{i}).body_speed.raw(o,1:length(tmp)) = speed4analysis(tmp);
                    elseif ix-(idx_Init-1) < 1 % first trials
                        tmp = 1:(ix + idx_End);
                        Stim_data.(Events{i}).idx(o,(n_trials-length(tmp)+1):end) = tmp;
                        for ii = 1:channel_number
                            Stim_data.(Events{i}).dFF.(dFF_names{ii}).raw(o,(n_trials-length(tmp)+1):end) = ...
                                streams.ZScoredFF.(dFF_names{ii})(tmp);
                            % first derivative:
                            dy = diff(Stim_data.(Events{i}).dFF.(dFF_names{ii}).raw(o,(n_trials-length(tmp)+1):end)); % differences between contiguous elements: substract point n to the point in n+1 --> 1 point smaller than original data
                            dx = diff(time_vect(tmp)); % differences between contiguous elements, time vector. doesnt matter if starts at 0. will have 1 element less                           
                            Stim_data.(Events{i}).firstderiv.(dFF_names{ii}).raw(o,(n_trials-length(tmp))+1:end-1) = dy./dx; % first derivative  - last point is a NaN
                        end
                        Stim_data.(Events{i}).body_speed.raw(o,(n_trials-length(tmp)+1):end) = speed4analysis(tmp);
                    end
                end
                for ii = 1:channel_number
                    Stim_data.(Events{i}).dFF.(dFF_names{ii}).baseline_corrected = ...
                        Stim_data.(Events{i}).dFF.(dFF_names{ii}).raw - nanmean...
                        (Stim_data.(Events{i}).dFF.(dFF_names{ii}).raw...
                        (:,t_trials >= BASELINE_PER(1) & ...
                        t_trials <= BASELINE_PER(2)),2);
                    Stim_data.(Events{i}).firstderiv.(dFF_names{ii}).baseline_corrected = ...
                        Stim_data.(Events{i}).firstderiv.(dFF_names{ii}).raw - nanmean...
                        (Stim_data.(Events{i}).firstderiv.(dFF_names{ii}).raw...
                        (:,t_trials >= BASELINE_PER(1) & ...
                        t_trials <= BASELINE_PER(2)),2); % added this (first derivative)
                end
                Stim_data.(Events{i}).body_speed.baseline_corrected = ...
                    Stim_data.(Events{i}).body_speed.raw - nanmean...
                    (Stim_data.(Events{i}).body_speed.raw...
                    (:,t_trials >= BASELINE_PER(1) & ...
                    t_trials <= BASELINE_PER(2)),2);
            end

            
            
            %% Calculations of dFF peaks, minima and latency thereof within 0-5sec window _ movement period, for each individual trial
               % initialize
            raw_or_corr ={'raw','baseline_corrected'};
            for r=1:length(raw_or_corr)
                for d=1:length(dFF_names) 
                    for i=1:length(Events)
                        Measurements.maxima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r}) = ones(size(Stim_data.(Events{i}).dFF.(dFF_names{d}).raw,1),1)*nan;
                        Measurements.minima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r}) = ones(size(Stim_data.(Events{i}).dFF.(dFF_names{d}).raw,1),1)*nan;
                        Measurements.minimalatency.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r}) = ones(size(Stim_data.(Events{i}).dFF.(dFF_names{d}).raw,1),1)*nan;
                        Measurements.maximalatency.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r}) = ones(size(Stim_data.(Events{i}).dFF.(dFF_names{d}).raw,1),1)*nan;
                        Measurements.maxima.(Events{i}).body_speed.(raw_or_corr{r}) = ones(size(Stim_data.(Events{i}).dFF.(dFF_names{d}).raw,1),1)*nan;
                        Measurements.minima.(Events{i}).body_speed.(raw_or_corr{r}) = ones(size(Stim_data.(Events{i}).dFF.(dFF_names{d}).raw,1),1)*nan;
                        Measurements.minimalatency.(Events{i}).body_speed.(raw_or_corr{r}) = ones(size(Stim_data.(Events{i}).dFF.(dFF_names{d}).raw,1),1)*nan;
                        Measurements.maximalatency.(Events{i}).body_speed.(raw_or_corr{r}) = ones(size(Stim_data.(Events{i}).dFF.(dFF_names{d}).raw,1),1)*nan;
                    end
                end
            end
            % allocate
            ix11 = find(t_trials==0);
            periodend = 1.5; % 1.5 seconds
            tp = find(t_trials-periodend>0);
            ix21 = tp(1);
            periodend = 2; % 2 seconds
            tp = find(t_trials-periodend>0);
            ix12 = tp(1);
            periodend = 4.5; % 4.5 seconds
            tp = find(t_trials-periodend>0);
            ix22 = tp(1);
            for r = 1:length(raw_or_corr)
                for d = 1:length(dFF_names)
                    for i = 1:length(Events)
                        for o = 1:size(Stim_data.(Events{i}).dFF.(dFF_names{d}).raw,1)
                            if i==1 % onset
                                Measurements.maxima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(o) = ...
                                    max(Stim_data.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(o,ix12:ix22),[],2);
                                Measurements.minima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(o) = ...
                                    min(Stim_data.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(o,ix11:ix21),[],2);
                                j=find(Stim_data.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(o,ix12:ix22) == Measurements.maxima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(o));
                                t_trials_temp=t_trials(ix12:ix22);
                                Measurements.maximalatency.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(o) = t_trials_temp(j);
                                k=find(Stim_data.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(o,ix11:ix21) == Measurements.minima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(o));
                                t_trials_temp=t_trials(ix11:ix21);
                                Measurements.minimalatency.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(o) = t_trials_temp(k);
                            else % offset
                                Measurements.maxima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(o) = ...
                                    max(Stim_data.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(o,ix11:ix21),[],2);
                                Measurements.minima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(o) = ...
                                    min(Stim_data.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(o,ix12:ix22),[],2);
                                j=find(Stim_data.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(o,ix11:ix21) == Measurements.maxima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(o));
                                t_trials_temp=t_trials(ix11:ix21);
                                if length(j) == 1
                                    Measurements.maximalatency.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(o) = t_trials_temp(j);
                                end
                                k=find(Stim_data.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(o,ix12:ix22) == Measurements.minima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(o));
                                t_trials_temp=t_trials(ix12:ix22);
                                if length(k) == 1
                                    Measurements.minimalatency.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(o) = t_trials_temp(k);
                                end
                            end
                        end      
                    end
                end
            end 
            for r = 1:length(raw_or_corr)
                for i = 1:length(Events)
                    for o = 1:size(Stim_data.(Events{i}).dFF.(dFF_names{d}).raw,1)
                        if i==1 % onset
                            Measurements.maxima.(Events{i}).body_speed.(raw_or_corr{r})(o) = ...
                                max(Stim_data.(Events{i}).body_speed.(raw_or_corr{r})(o,ix12:ix22),[],2);
                            Measurements.minima.(Events{i}).body_speed.(raw_or_corr{r})(o) = ...
                                min(Stim_data.(Events{i}).body_speed.(raw_or_corr{r})(o,ix11:ix21),[],2);
                            j=find(Stim_data.(Events{i}).body_speed.(raw_or_corr{r})(o,ix12:ix22) == Measurements.maxima.(Events{i}).body_speed.(raw_or_corr{r})(o));
                            t_trials_temp=t_trials(ix12:ix22);
                            Measurements.maximalatency.(Events{i}).body_speed.(raw_or_corr{r})(o) = t_trials_temp(j);
                            k=find(Stim_data.(Events{i}).body_speed.(raw_or_corr{r})(o,ix11:ix21) == Measurements.minima.(Events{i}).body_speed.(raw_or_corr{r})(o));
                            t_trials_temp=t_trials(ix11:ix21);
                            Measurements.minimalatency.(Events{i}).body_speed.(raw_or_corr{r})(o) = t_trials_temp(k);
                        else % offset
                            Measurements.maxima.(Events{i}).body_speed.(raw_or_corr{r})(o) = ...
                                max(Stim_data.(Events{i}).body_speed.(raw_or_corr{r})(o,ix11:ix21),[],2);
                            Measurements.minima.(Events{i}).body_speed.(raw_or_corr{r})(o) = ...
                                min(Stim_data.(Events{i}).body_speed.(raw_or_corr{r})(o,ix12:ix22),[],2);
                            j=find(Stim_data.(Events{i}).body_speed.(raw_or_corr{r})(o,ix11:ix21) == Measurements.maxima.(Events{i}).body_speed.(raw_or_corr{r})(o));
                            t_trials_temp=t_trials(ix11:ix21);
                            if length(j) == 1
                                Measurements.maximalatency.(Events{i}).body_speed.(raw_or_corr{r})(o) = t_trials_temp(j);
                            end
                            k=find(Stim_data.(Events{i}).body_speed.raw(o,ix12:ix22) == Measurements.minima.(Events{i}).body_speed.(raw_or_corr{r})(o));
                            t_trials_temp=t_trials(ix12:ix22);
                            if length(k) == 1
                                Measurements.minimalatency.(Events{i}).body_speed.(raw_or_corr{r})(o) = t_trials_temp(k);
                            end
                        end
                    end      
                end
            end
       
            % plot
            for r=1:length(raw_or_corr)
                for d=1:length(dFF_names)
                    figure; 
                    for i=1:length(Events)
                        if i == 1; u = 1; else; u = 3; end;
                        subplot(2,2,u)
                        plot(Measurements.minimalatency.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r}),'Linewidth',2);
                        hold on;
                        plot(Measurements.minimalatency.(Events{i}).body_speed.(raw_or_corr{r}),'Linewidth',2);
                        title([dFF_names{d},', mov ',Events{i}(5:end),' latency to min']);
                        legend('dFF latency','speed latency');
                        ylabel('latency (s)')
                        if i == 1; u = 2; else; u = 4; end; 
                        subplot(2,2,u)
                        plot(Measurements.maximalatency.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r}),'Linewidth',2);
                        hold on;
                        plot(Measurements.maximalatency.(Events{i}).body_speed.(raw_or_corr{r}),'Linewidth',2);
                        title([dFF_names{d},', mov ',Events{i}(5:end),' latency to max']);
                        legend('dFF latency','speed latency');
                        ylabel('latency (s)')
                        if r == 1; sgtitle('raw'); else; sgtitle('baseline corrected'); end
                    end
                    if save_plot
                        saveas(gcf,[PATH2SAVE,'figures\Measurements latency',dFF_names{d},' ',raw_or_corr{r}],'tif');
                        saveas(gcf,[PATH2SAVE,'figures\Measurements latency',dFF_names{d},' ',raw_or_corr{r}],'fig');
                    end
                end
            end            
                
            

            %% Calculations of dFF peaks, minima and latency thereof within 0-5sec window: baseline epoch
            % initialize
            for r=1:length(raw_or_corr)
                for d=1:length(dFF_names) 
                    for i=1:length(Events)
                        Measurements_con.maxima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r}) = ones(size(Stim_data.(Events{i}).dFF.(dFF_names{ii}).raw,1),1)*nan;
                        Measurements_con.minima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r}) = ones(size(Stim_data.(Events{i}).dFF.(dFF_names{ii}).raw,1),1)*nan;
                        Measurements_con.maxima.(Events{i}).body_speed.(raw_or_corr{r}) = ones(size(Stim_data.(Events{i}).dFF.(dFF_names{ii}).raw,1),1)*nan;
                        Measurements_con.minima.(Events{i}).body_speed.(raw_or_corr{r}) = ones(size(Stim_data.(Events{i}).dFF.(dFF_names{ii}).raw,1),1)*nan;
                    end
                end
            end
            % allocate
            periodend = -4.5; % 1.5 seconds
            tp = find(t_trials-periodend>0);
            ix1 = tp(1);
            periodend = -2.5; % 2 seconds
            tp = find(t_trials-periodend>0);
            ix2 = tp(1);
            for r = 1:length(raw_or_corr)
                for d = 1:length(dFF_names)
                    for o = 1:size(Stim_data.(Events{i}).dFF.(dFF_names{d}).raw,1)
                        i=1; % onset
                        Measurements_con.maxima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(o) = ...
                            max(Stim_data.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(o,ix1:ix2),[],2);
                        i=2; % offset
                        Measurements_con.minima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(o) = ...
                            min(Stim_data.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(o,ix1:ix2),[],2);
                    end      
                end
            end 
            for r = 1:length(raw_or_corr)
                for o = 1:size(Stim_data.(Events{i}).dFF.(dFF_names{d}).raw,1)
                    i=1; 
                    Measurements_con.maxima.(Events{i}).body_speed.(raw_or_corr{r})(o) = ...
                        max(Stim_data.(Events{i}).body_speed.(raw_or_corr{r})(o,ix1:ix2),[],2);
                    i=2;
                    Measurements_con.minima.(Events{i}).body_speed.(raw_or_corr{r})(o) = ...
                        min(Stim_data.(Events{i}).body_speed.(raw_or_corr{r})(o,ix1:ix2),[],2);
                end      
            end
       
            % plot the minima and maxima of control and baseline
            for r=1:length(raw_or_corr)
                figure; 
                for d=1:length(dFF_names)
                        i=1; % onset, look at maxima
                        if d == 1; u = 1; else; u = 3; end;
                        subplot(2,2,u)
                        plot(Measurements.maxima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r}),'Linewidth',2);
                        hold on;
                        plot(Measurements_con.maxima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r}),'Linewidth',2);                        
%                         plot(Measurements.minima.(Events{i}).body_speed.(raw_or_corr{r}),'Linewidth',2);
                        title([dFF_names{d},', mov ',Events{i}(5:end),' maxima']);
                        legend('movement','control');
                        ylabel('dFF maxima')
                        i=2; % offset, look at maxima
                        if d == 1; u = 2; else; u = 4; end;                        
                        subplot(2,2,u)
                        plot(Measurements.minima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r}),'Linewidth',2);
                        hold on;
                        plot(Measurements_con.minima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r}),'Linewidth',2);                       
%                         plot(Measurements.maximalatency.(Events{i}).body_speed.(raw_or_corr{r}),'Linewidth',2);
                        title([dFF_names{d},', mov ',Events{i}(5:end),' minima']);
                        legend('movement','control');
                        ylabel('dFF minima')
                        if r == 1; sgtitle('raw'); else; sgtitle('baseline corrected'); end
                    if save_plot
                        saveas(gcf,[PATH2SAVE,'figures\Measurements vs baseline',dFF_names{d},' ',raw_or_corr{r}],'tif');
                        saveas(gcf,[PATH2SAVE,'figures\Measurements vs baseline',dFF_names{d},' ',raw_or_corr{r}],'fig');
                    end
                end
            end            
             
            
            %% plot average
             % plot the minima and maxima of control and baseline
            color2plot_new = {[0.72,0.27,1.00],[0.52,0.27,1.00],[0.47,0.9,0.19],[0,0.7,0.19]}; % purple, purple-ish, green, green-ish

            for r=1:length(raw_or_corr)
                figure; 
                i=1; % onset, look at maxima
                subplot(1,2,i)
                for d=1:length(dFF_names)
                    %control
                    if d == 1; z = 1; else; z = 3; end
                    tmp_avg = nanmean(Measurements_con.maxima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r}),1);
                    tmp_error = nanstd(Measurements_con.maxima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r}),1,1)./...
                        sqrt(sum(~isnan(Measurements_con.maxima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(:,1))));
%                     error_area([],tmp_avg,tmp_error,color2plot{1},0.25); % error_area(X, Y, barlength, color, alpha, varargin) to plot SEM  
                    bar([z],tmp_avg,'FaceColor',color2plot_new{z},'EdgeColor',color2plot_new{z},'LineWidth',1.5);
                    hold on;
                    er = errorbar([z],tmp_avg,tmp_error,'-'); % errorbar(x,y,err) % color2plot{1},0.25
                    er.Color = [0 0 0]; er.LineStyle = 'none'; er.LineWidth = 2;
                    %movement
                    if d == 1; z = 2; else; z = 4; end                    
                    tmp_avg2 = nanmean(Measurements.maxima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r}),1);
                    tmp_error2 = nanstd(Measurements.maxima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r}),1,1)./...
                        sqrt(sum(~isnan(Measurements.maxima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(:,1))));
%                     error_area([],tmp_avg,tmp_error,color2plot{1},0.25); % error_area(X, Y, barlength, color, alpha, varargin) to plot SEM  
                    bar([z],tmp_avg2,'FaceColor',color2plot_new{z},'EdgeColor',color2plot_new{z},'LineWidth',1.5);
                    hold on;
                    er = errorbar([z],tmp_avg2,tmp_error2,'-'); % errorbar(x,y,err) % color2plot{1},0.25
                    er.Color = [0 0 0]; er.LineStyle = 'none'; er.LineWidth = 2;
                end
                title(['mov ',Events{i}(5:end),', maxima']);
                legend('GPe con',' ','GPe mvt',' ','SNr con',' ','SNr mvt',' ','Location','Northeast','NumColumns' , 2);
                ylabel('dFF maxima')
                min_x=0; max_x=5; min_y=0; max_y=2; axis([min_x max_x min_y max_y]); 
                set(gca,'xtick',[])

                i=2; % offset, look at minima
                subplot(1,2,i)
                for d=1:length(dFF_names)
                    %control
                    if d == 1; z = 1; else; z = 3; end
                    tmp_avg = nanmean(Measurements_con.minima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r}),1);
                    tmp_error = nanstd(Measurements_con.minima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r}),1,1)./...
                        sqrt(sum(~isnan(Measurements_con.minima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(:,1))));
%                     error_area([],tmp_avg,tmp_error,color2plot{1},0.25); % error_area(X, Y, barlength, color, alpha, varargin) to plot SEM  
                    bar([z],tmp_avg,'FaceColor',color2plot_new{z},'EdgeColor',color2plot_new{z},'LineWidth',1.5);
                    hold on;
                    er = errorbar([z],tmp_avg,tmp_error,'-'); % errorbar(x,y,err) % color2plot{1},0.25
                    er.Color = [0 0 0]; er.LineStyle = 'none'; er.LineWidth = 2;
                    %movement
                    if d == 1; z = 2; else; z = 4; end                    
                    tmp_avg2 = nanmean(Measurements.minima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r}),1);
                    tmp_error2 = nanstd(Measurements.minima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r}),1,1)./...
                        sqrt(sum(~isnan(Measurements.minima.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(:,1))));
%                     error_area([],tmp_avg,tmp_error,color2plot{1},0.25); % error_area(X, Y, barlength, color, alpha, varargin) to plot SEM  
                    bar([z],tmp_avg2,'FaceColor',color2plot_new{z},'EdgeColor',color2plot_new{z},'LineWidth',1.5);
                    hold on;
                    er = errorbar([z],tmp_avg2,tmp_error2,'-'); % errorbar(x,y,err) % color2plot{1},0.25
                    er.Color = [0 0 0]; er.LineStyle = 'none'; er.LineWidth = 2;
                end
                title(['mov ',Events{i}(5:end),', minima']);
                legend('GPe con',' ','GPe mvt',' ','SNr con',' ','SNr mvt',' ','Location','SouthEast','NumColumns' , 2);
                ylabel('dFF minima')
                min_x=0; max_x=5; min_y=-1; max_y=0; axis([min_x max_x min_y max_y]); 
                set(gca,'xtick',[])

                if r == 1; sgtitle('raw'); else; sgtitle('baseline corrected'); end
                if save_plot
                    saveas(gcf,[PATH2SAVE,'figures\Measurements vs baseline ave',dFF_names{d},' ',raw_or_corr{r}],'tif');
                    saveas(gcf,[PATH2SAVE,'figures\Measurements vs baseline ave',dFF_names{d},' ',raw_or_corr{r}],'fig');
                end
            end         
            
            
            %% Average plot of speed only
            raw_or_corr = {'raw','baseline_corrected'};
%             color2plot = {'b','g','m','r'};
            color2plot = {[1.00,0.41,0.16],[0.72,0.27,1.00],[0.47,0.9,0.19],[1.0 0.1 0.3],[0 0 0.7],[0.5 0.5 0.5]}; % orange, purple, green,red, blue, grey 

            % speed
            for r=1:length(raw_or_corr)
%             for r=2
                for i=1:length(Events)                  
                        figure; 
                        if show_plot == 0
                            set(gcf,'visible','off')
                        end
                        hold on
                        tmp_avg = nanmean(Stim_data.(Events{i}).body_speed.(raw_or_corr{r}),1);
                        tmp_error = nanstd(Stim_data.(Events{i}).body_speed.(raw_or_corr{r}),1,1)./...
                            sqrt(sum(~isnan(Stim_data.(Events{i}).body_speed.(raw_or_corr{r})(:,1))));
                        error_area(t_trials,tmp_avg,tmp_error,color2plot{1},0.25); % error_area(X, Y, barlength, color, alpha, varargin) to plot SEM
%                         max_avg(pow) = max(tmp_avg);
                        xline(0,'-k'); %
                        xline(stim_duration,'-m');
                        yline(0,'-.k');
%                         yline(max(max_avg),'-.r');
                        xlabel('Time (s)');
                        ylabel('Body speed (cm)');
%                         if i == 4; % for speed
%                             ylabel([datatype{i},' m/s']);
%                         end
                        xlim([TRANGE(1) TRANGE(2)]);
%                         ylim(limits2plot.(raw_or_corr{r}).(dFF_names{d}).(datatype{i})); %see parameters setting at start                     
%                         PeakAve = num2str(mean(max_avg));
%                         PeakMax = num2str(max(max_avg));
                        sgtitle([raw_or_corr{r},' ',Events{i},' ','speed'],'Interpreter','none')
%                         sgtitle([raw_or_corr{r},' ',datatype{i},' ',dFF_names{d},', AvePeak: ',PeakAve,', MaxPeak: ',PeakMax],'Interpreter','none')

                        % legend
%                         legend([Opto_Powers{1}(2:end),'error'],[Opto_Powers{1}(2:end)],...
%                                 [Opto_Powers{2}(2:end),'error'],[Opto_Powers{2}(2:end)],...
%                                 'Opto On','Opto off','Location','northwest','NumColumns',1)      

                        %saveplot or not
                        if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                            saveas(gcf,[PATH2SAVE,'figures\AVE movement speed ',raw_or_corr{r},'.tif'])
                            saveas(gcf,[PATH2SAVE,'figures\AVE movement speed ',raw_or_corr{r},'.fig'])
                        end
                    end
            end
                
                
            
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
                        tmp_avg = nanmean(Stim_data.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r}),1);
                        tmp_error = nanstd(Stim_data.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r}),1,1)./...
                            sqrt(sum(~isnan(Stim_data.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(:,1))));
                        error_area(t_trials,tmp_avg,tmp_error,color2plot{2},0.25); % error_area(X, Y, barlength, color, alpha, varargin) to plot SEM
                        % dFF SNr
                        d=2;
%                         yyaxis left
                        tmp_avg = nanmean(Stim_data.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r}),1);
                        tmp_error = nanstd(Stim_data.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r}),1,1)./...
                            sqrt(sum(~isnan(Stim_data.(Events{i}).dFF.(dFF_names{d}).(raw_or_corr{r})(:,1))));
                        error_area(t_trials,tmp_avg,tmp_error,color2plot{3},0.25); % error_area(X, Y, barlength, color, alpha, varargin) to plot SEM

                           if i==1
                            if r==2
                                min_x=0; max_x=5; min_y=-0.5; max_y=1.5; axis([min_x max_x min_y max_y]); 
                            else
                                min_x=0; max_x=5; min_y=-0.8; max_y=1; axis([min_x max_x min_y max_y]);                                 
                            end
                        else
                            if r==2
                                min_x=0; max_x=5; min_y=-1.5; max_y=0.5; axis([min_x max_x min_y max_y]);                             
                            else
                                min_x=0; max_x=5; min_y=-0.8; max_y=1.2; axis([min_x max_x min_y max_y]);                                                             
                            end
                        end
                        % speed
                        yyaxis right
                        tmp_avg3 = nanmean(Stim_data.(Events{i}).body_speed.(raw_or_corr{r}),1);
                        tmp_error3 = nanstd(Stim_data.(Events{i}).body_speed.(raw_or_corr{r}),1,1)./...
                            sqrt(sum(~isnan(Stim_data.(Events{i}).body_speed.(raw_or_corr{r})(:,1))));
                        error_area(t_trials,tmp_avg3,tmp_error3,color2plot{1},0.25);
                        ylabel('Mouse speed (cm/s)'); 
                        if i==1
                            if r==2
                                min_x=0; max_x=5; min_y=-2; max_y=8; axis([min_x max_x min_y max_y]); 
                            else
                                min_x=0; max_x=5; min_y=1; max_y=8; axis([min_x max_x min_y max_y]); 
                            end
                        else
                           if r==2
                                min_x=0; max_x=5; min_y=-8; max_y=2; axis([min_x max_x min_y max_y]); 
                           else
                                min_x=0; max_x=5; min_y=0; max_y=10; axis([min_x max_x min_y max_y]);                                
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
%                               sgtitle([raw_or_corr{r}],'Interpreter','none')
%                         sgtitle([raw_or_corr{r},' ',datatype{i},' ',dFF_names{d},', AvePeak: ',PeakAve,', MaxPeak: ',PeakMax],'Interpreter','none')

                        % legend
                        legend('dFF GPe',' ','dFF SNr',' ','Mouse speed',' ','Location','Southwest','NumColumns',1)      
                        legend boxoff
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
            
          

        end % end of if loop "if done"
 
          
        
        
        
        %% SAVE INDIVDIDUAL SESSION   
        
        if BehaviorDeeplabcut == 0
            resample_Position = {};
        end
        
        if done == 0 || overwrite == 1
            Params.MICEId = mice_list(nummice).name;
            Params.SESSIONId = sessions(s).name;
            Params.Virus = Virus;
            Params.ExperimentType = ExperimentType;
            Params.Preprocessing.lowpassfiltering = lowpassfiltering;
            Params.Preprocessing.lowpassfreq = lowpassfreq;
            Params.Preprocessing.detrending = detrending;
            Params.Preprocessing.deletemovingwindow = deletemovingwindow;
            Params.Preprocessing.highpassfiltering = highpassfiltering;
            Params.Preprocessing.highpassfreq = highpassfreq;
            Params.Path.PATH2SESSION = PATH2SESSION;
            Params.Path.PATH2SAVE = PATH2SAVE;          
            save([PATH2SAVE,'IndividualData.mat'],'Params','body_mov','Stim_data','allcorr','resample_Position','epocs','datatype','dFF_names',...
                'dt_ds','sampling_rate_ds','length_data','streams','time_vect','TRANGE','BASELINE_WIN','t_trials','n_trials','Events','BodySpeed','raw_or_corr',...
                'Measurements','Measurements_con'); 
        end    


        %% Ending all the loops

    end %for s=1:length sessions
end % for all mice: i=1:numfiles, see above

    
