% DATA ANALYSIS PIPELINE FOR FIBER PHOTOMETRY EXPERIMENT: OPEN FIELD WITH OPTO STIM: 
% Marie_FP_IndivData_extraction 
% Closed loop experiment based on dFF threshold or motor bout. 
% also calculates average speed and average rotations (in the stim epoch

% extracts FP data and behavioral data 
% generates and saves a matrix with individual trials aligned to specific event
%       for now events are optogenetic stimulation
% generates and saves graphs aligned to specific event


%% INITIALIZATIONS
close all; clear all; clc;
set(0,'defaultfigurecolor',[1 1 1])

%% SETUP PARAMETERS
%%%%%%%%%% TTL EXTRACTION
ttl_extract = 'default'; %'default' or 'special' (needs to be 7 characters)
 
%%%%%%%%%% DATA EXTRACTION
channel_number = 1; % 1 if 1-site recording, 2 if 2-site recording
channel_names = {'GPe'}; % eg brain region
color_number = 2; % 2 if 405 and 465 recording, 3 if also 565
color_names = {'c405','c465'}; % can also be 565
dFF_number = 1; % how many dFF to analyze, can be 1 or 2, based on number of channels and colors (if more than 2, need to edit script)
dFF_names = {'GPe'}; %put them in this order: 1) channel 1, 465 and 2) channel 1, 560 OR channel 2: 465 (if other combinations, need to edit script)
datatype = {'ZScoredFF','dFF','dFFsy'} %,'speed'}; %data analyses of interest -- remove speed if not analysing behavior
% datatype = {'speed','angle','cumangle'}; %; %data analyses of interest -- remove speed if not analysing behavior

calcdFF = 1; % 1 if you dont want to calculate dFF, 0 otherwise

%%%%%%%%%%% PREPROCESSING
% trimming
timetrim_start_set = 60; %time to chop off at start of recording, in seconds: adjust this manually, eg 60 seconds
timetrim_end_set = 1; %time to chop off at end of recording, eg 1 second
% low pass filtering of raw 405 and 470
lowpassfiltering = 0; % 1 if you want to low pass filter- typically yes- 0 for no
lowpassfreq = 1; % typically 1Hz
% detrending dFF
detrending = 0;% %% Write 1 for detrend dFF, normal detrend function of Matlab- typically no (Write 0) for open field
% substract 8th percentile as a baseline calculate using a moving window
deletemovingwindow = 0; % 1 to delete moving baseline 8th percentile, 0 if no. Typically no for open field, yes for rotarod
% high pass filtering of dFF
highpassfiltering = 0; %'yes' if you want to high pass filter, 0 for no. Typically no for open field
highpassfreq = 0.005; % typically 0.005Hz

%%%%%%%%%%% TRIAL DEFINITION
% time window for the trials and graphs
TRANGE = [-15 15]; % will create events for a -5 to +5 sec window (you can always plot less later)
% baseline correction window
BASELINE_WIN.OptoStim = [-15 5];% baseline correction window% BASELINE_PER = [-5 -1]; 
% variables to align to
Events2Align2 = {'PulseStart_time','PulseStop_time'}; %can be {'optostim_onset','optostim_offset'} for opto, or any other TTL epoc in the "epocs" structure with same organization

%%%%%%%%%%% EXPERIMENT TYPE
ExperimentType = 'power'; % Options should be 5 characters: 'power', 'freq_', 'mobil','poLED' (for photometryLED). If other, then the trials are not defined and need to go in and change trial allocation.
BehaviorAnymaze = 0; % 0 if not, 1 if FP rig received frames from Anymaze interface and stored its data into a CSV file (1 row per frame)
BehaviorDeeplabcut = 1; % 0 if not, 1 if FP took camera frames itself and it will be analyzed by another program and results put into a CSV file (1 row per frame)
cumanglecalc = 0; % calculate cumulative angle or not
AnalyzeDeepLabCut = 0; % analyze DLC open field data or not

Opto_Powers = {'p0uW','p500uW','p2000uW'}; % POWER DEPENDENCE, dFF
anymazepower = 0.2; % 0,2mW for the mobility experiments
TrialType_Number = length(Opto_Powers);
stim_duration = 10; % 10 seconds or 3 seconds or 5 seconds
stim_duration2 = [];  % in case their is a ramp or another duration of interest

%%%%%%%%%%% HOW MANY MICE TO ANALYZE
LoopOrNot = 1; % 1 to loop, 0 if you only want to test one mouse and dont want to loop --> if so, edit the mouse number you want below: nummice_set
nummice_set = 4;

%%%%%%%%%%% SHOW, SAVE, OVERWRITE PARAMETERS
set(gcf,'visible','off')
show_plot = 1; % If 0, plots are not displayed 
save_plot = 0; % If 0, plots are not saved
reanalysis = 1; % If 1, the code runs for sessions already analyzed. If 0, session excluded if analysed. This is based on the existence of the matlab space "IndividualData.mat"
overwrite = 1; % If reanalyzing data, 1 if you want to save them and overwrite results and 0 if not          
Collapse_sessions = 0; % If 0, analysis only of individual sessions. If 1, analysis of the group data (pooled)          
Color_scale = []; % For the heatmaps. If empty, automatically adjusted to the data.   


%% DEFINE PATHS
path2data = uigetdir('select folder'); % SELECT FOLDER OF GROUP TO ANALYZE, eg Chrimson or mCherry folder - %Location of the data
mice_list = dir(path2data); %all things in this folder
% Define paths and data to analyze
path2savefolder = path2data; %Path to save


%% IDENTIFY MICE TO ANALYZE 
for o = length(mice_list):-1:1
    if mice_list(o).isdir == 0  %remove non-folders
        mice_list(o) = [];
    else
        if  strcmp(mice_list(o).name,'data') == 1 || strcmp(mice_list(o).name,'figures') == 1 ...   %remove folders with data or figures or ., ..
            || contains(mice_list(o).name,'data') || contains(mice_list(o).name,'figures') || contains(mice_list(o).name,'other') ...
            || contains(mice_list(o).name,'results') || contains(mice_list(o).name,'turns') || strcmp(mice_list(o).name,'.') == 1 || strcmp(mice_list(o).name,'..') == 1;
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
    Virus = Virus_cell{:}; %just extracting the string out of the cell array   --> "Chrimson or mCherry (name of the folder you selected)
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
        elseif strcmp(sessions(o).name,'data') == 1 || strcmp(sessions(o).name,'figures') == 1 || strcmp(sessions(o).name,'results') == 1 ...
             || strcmp(sessions(o).name,'.') == 1 || strcmp(sessions(o).name,'..') == 1 || contains(sessions(o).name,'.csv') == 1 ...
             || contains(sessions(o).name,'data') == 1 || contains(sessions(o).name,'figures') == 1 || contains(sessions(o).name,'other') == 1 || contains(sessions(o).name,'BACKUP') ...
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
            stream_name_type405A = {'x05A','05A','A05A'}; %possible names on kellendonk/ansorge rigs- can add to this list
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
            sampling_rate = data.streams.(stream_name_type405A{1}).fs;
            N = 10; %downsample 10 times (from initial sampling rate of 1017.3 Hz) 
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

            
            %% ANYMAZE: %MATLAB SLOW ON THIS
            % Transforming camera TTLs/frames that are beyond the end of the time vector ,camera takes longer to turn off.
            if BehaviorAnymaze == 1 
                % last time point of recording in seconds
                time_vect_end = time_vect(length(time_vect)); 
                    frontcam_timevect_lessthan = epocs.Framesync_start_time(epocs.Framesync_start_time <= time_vect_end); 
                FrontCamOnLocations = zeros(1,length(frontcam_timevect_lessthan)); 
                for i = 1:length(frontcam_timevect_lessthan)
                    if epocs.Framesync_start_time(i) <= time_vect_end 
                        tmp = find(time_vect >= epocs.Framesync_start_time(i));  
                        FrontCamOnLocations(i) = tmp(1); 
                    end
                end

                epocs.FrontCamTTL = zeros(1,length_data);
                for i=1:length(FrontCamOnLocations)
                    epocs.FrontCamTTL(FrontCamOnLocations(i)-1) = 5; 
                end
                clear frontcam_timevect_lessthan
                
                %% LOAD THE BEHAVIORAL INFORMATION
                % Here this is formatted for uploading Anymaze data in a CSV file, but can easily be adjusted for other CSV
                % behavioral data, just need to adjust according to CSV format as well as TTL alignment below
                files_csv = dir(fullfile(PATH2SESSION,'*.csv'));
                numfiles = length(files_csv); 
                tmp_name = files_csv.name(1:end);
                tmp_table = readtable([files_csv.folder,'\',tmp_name]); %readtable: create table from file which is a csv file in this folder with the same name
                Body_parts = tmp_table.Properties.VariableNames(2:end); %body parts can be the position of a specific body part, or can be a behavioral variable like speed

                for k = 1:size(Body_parts,2)
                    % Behavioral data
                    tmp_array = table2array(tmp_table(1:end,k+1));
                    BehavData.(AnimalID).(Body_parts{k}) = tmp_array;
                end


                %% CHECK CONSISTENCY BETWEEN CAM TTLs AND VIDEO DATA
                if ~isempty(Body_parts)
                    FP_TTL_number = length(find(epocs.FrontCamTTL)) %leave visible so I can compare
                    Behav_frame_number = length(BehavData.(AnimalID).Speed) %leave visible so I can compare
                        if Behav_frame_number > FP_TTL_number
                        for i = 1:length(Body_parts)
                            BehavData.(AnimalID).(Body_parts{i})(FP_TTL_number+1:end) = [];
                        end
                    end
                    if Behav_frame_number < FP_TTL_number
                        for i = 1:length(Body_parts)
                            BehavData.(AnimalID).(Body_parts{i})(FP_TTL_number) = NaN;
                            BehavData.(AnimalID).(Body_parts{i})(FP_TTL_number+1:end) = NaN;

                        end
                    end
                end

                FP_TTL_number_post = length(find(epocs.FrontCamTTL)) %leave visible so I can compare
                Behav_frame_number_post = length(BehavData.(AnimalID).Speed) %leave visible so I can compare

                %% CREATE NEW STRUCTURE NewPos WITH BODY POSITION THE SAME SIZE AS FP DATA
                if ~isempty(Body_parts)
                    for i = 1:length(Body_parts)
                        newPos.(AnimalID).(Body_parts{i}) = nan(1,length_data); %create new vector of nans of the size of the body parts (Eg 5 of them) and the size of data_fp (number of frames)
                    end
                    for i = 1:length(Body_parts)
                        newPos.(AnimalID).(Body_parts{i})(epocs.FrontCamTTL(1,:) == 5) = BehavData.(AnimalID).(Body_parts{i}); %put into newPos the behav data. %but only every time the TTL is equal to 5 (TTL);
                    end
                end
            end %end of Anymaze/Behavioral data extraction


            
            
            
            %% DEEPLABCUT
            if BehaviorDeeplabcut == 1 
                if AnalyzeDeepLabCut == 1
                % last time point of recording in seconds
                time_vect_end = time_vect(length(time_vect)); 
                     frontcam_timevect_lessthan = epocs.FrontCam_start_time(epocs.FrontCam_start_time <= time_vect_end); 
                FrontCamOnLocations = zeros(1,length(frontcam_timevect_lessthan)); 
                for i = 1:length(frontcam_timevect_lessthan)
                    if epocs.FrontCam_start_time(i) <= time_vect_end 
                        tmp = find(time_vect >= epocs.FrontCam_start_time(i));  
                        FrontCamOnLocations(i) = tmp(1); 
                    end
                end

                epocs.FrontCamTTL = zeros(1,length_data);
                for i=1:length(FrontCamOnLocations)
                    epocs.FrontCamTTL(FrontCamOnLocations(i)-1) = 5;    
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
                set(gcf,'visible','off')

                f= figure; f.Position = [50 70 1800 900]; sgtitle('Likelihood for each body part across time before correction')
                for cam = 1:length(cameras) %:length(cameras)
                    tmpfld = fieldnames(Position{cam});
                    for i = 1:length(tmpfld)
                        idx2extract = i;      
                        subplot(2,round(length(tmpfld)./2),i); plot(Position{cam}.(tmpfld{idx2extract}).likelihood)
                        title([tmpfld{idx2extract}])
                    end
                end
%                 set(gcf,'visible','on')

                %interpolation
                set(gcf,'visible','off')
               
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
%                 set(gcf,'visible','on')
  
                %likelihood
                set(gcf,'visible','off')

                f= figure; f.Position = [50 70 1800 900]; sgtitle('Likelihood for each body part across time after interpolation correction (keep point if likelihood >0.9)')
                for cam = 1:length(cameras) %:length(cameras)
                    tmpfld = fieldnames(Position{cam});
                    for i = 1:length(tmpfld)
                        idx2extract = i;      
                        subplot(2,round(length(tmpfld)./2),i); plot(Position{cam}.(tmpfld{idx2extract}).likelihood)
                        title([tmpfld{idx2extract}])
                    end
                end
%                 set(gcf,'visible','on')

                            
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
                   
            end %end of DLC Behavioral data extraction
            
       
            end
            
            
            %% Clear TDT data structure (now everything we need is in "streams" or "epocs")
%             clear data

            
            %% TRIMMING
            % Setting up trim indexes
            if BehaviorDeeplabcut == 1 
                timetrim_start = timetrim_start_set*1; 
                timetrim_end =  timetrim_end_set*1; 
                dummie1 = 1:length(time_vect); 
                dummie2 = dummie1(time_vect > timetrim_start); %only keep the indexes starting from the new start
                idx_start = dummie2(1); %index for the new start
                dummie2 = dummie1(time_vect > time_vect(end) - timetrim_end); %only keep the indexes from the new end to the current end
                idx_end = dummie2(1); %index for the new end
                clear dummie1 dummie2 
                
            elseif BehaviorAnymaze == 1 
                % trim end
                FrontCamTTL_on=find(epocs.FrontCamTTL); % indexes when FrontCamTTL not 0
                timetrim_end = time_vect(end)-((FrontCamTTL_on(end)+1)*dt_ds); % 
                dummie1 = 1:length(time_vect);
                dummie2 = dummie1(time_vect > time_vect(end) - timetrim_end);
                idx_end = dummie2(1); %to find the exact index within time_vect
                timetrim_end = time_vect(idx_end); %this is the time that corresponds to the actual index where you trim off
                timetrim_start = timetrim_start_set*1; %set at the beginning of the script- first delete eg first 90 seconds
                dummie2 = dummie1(time_vect > timetrim_start);
                idx_start = dummie2(1); %to find the exact index within time_vect- temporary; then will correct below, based on camera frames
                if ~isempty(Body_parts)
                    FrontCamTTL_trim = epocs.FrontCamTTL(idx_start:end);
                    FrontCamTTL_trim_on = find(FrontCamTTL_trim); 
                    idx_start = idx_start + FrontCamTTL_trim_on(1)-1; 
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
            if BehaviorAnymaze == 1
                if ~isempty(Body_parts)
                    for i = 1:length(Body_parts)
                        newPos_trim.(AnimalID).(Body_parts{i}) = newPos.(AnimalID).(Body_parts{i})(:,idx_start:idx_end); %put into newPos the behav data. %but only every time the TTL is equal to 5 (TTL);
                    end
                end
            end
            if BehaviorDeeplabcut == 1
                if AnalyzeDeepLabCut == 1
                for i = 1:length(tmpfld)
                    for ii = 1:length(tmpfld2)
                        newPos_trim{cam}.(tmpfld{i}).(tmpfld2{ii}) = newPos{cam}.(tmpfld{i}).(tmpfld2{ii})(:,idx_start:idx_end);
                    end
                end
                end
            end
            
            % Trimming: Plot
            for channel=1:length(channel_names)
                for colore=1:length(color_names)
                    figure; clf; 
%                     if show_plot == 0
                       set(gcf,'visible','off')
%                     end
                    %Plot after
                    subplot(2,1,1);
                    plot(time_vect,streams.rawdata.(channel_names{channel}).(color_names{colore})); hold on; 
                    for o = 1:size(epocs.PulseStart_time)
                        plot([epocs.PulseStart_time(o) epocs.PulseStart_time(o)],[-4 12],'m')
                        plot([epocs.PulseStop_time(o) epocs.PulseStop_time(o)],[-4 12],'c')
                    end
                    min_x=-200; max_x=max(time_vect)+200; min_y=mean(streams.rawdata.(channel_names{channel}).(color_names{colore}))-100; 
                    max_y=mean(streams.rawdata.(channel_names{channel}).(color_names{colore}))+100; axis([min_x max_x min_y max_y]); 
                    xlabel('time(sec)'); ylabel('fluorescence'); title('Data before trimming'); 
                    xline(timetrim_start,':k','trim'); xline(time_vect(end)-timetrim_end,':k','trim'); 
                    %Plot after
                    subplot(2,1,2); plot(time_vect_trim,streams.rawdata_trim.(channel_names{channel}).(color_names{colore})); hold on; 
                    for o = 1:size(epocs_trim.PulseStart_time)
                        plot([epocs_trim.PulseStart_time(o) epocs_trim.PulseStart_time(o)],[-4 12],'m')
                        plot([epocs_trim.PulseStop_time(o) epocs_trim.PulseStop_time(o)],[-4 12],'c')
                    end
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
            if BehaviorAnymaze == 1
                if ~isempty(Body_parts)
                    BehavData = newPos_trim;
                end
            end
            if BehaviorDeeplabcut == 1
                if AnalyzeDeepLabCut == 1
                for cam = 1:length(cameras)
                    if ~isempty(newPos_trim{cam})
                        % Ytail = Ytail_trim;
                        Position{cam} = newPos_trim{cam};
                    end
                end
                clear newPos_trim
                end
            end
            length_data = length(streams.rawdata.(channel_names{channel}).(color_names{colore}));
            clear epocs_trim time_vect_trim newPos_trim new_Pos
            streams = rmfield(streams,'rawdata_trim'); %clear variable within structure.

            if AnalyzeDeepLabCut == 1
            % Need to do FrontCamTTL again after trimming
                time_vect_end = time_vect(length(time_vect)); 
                    frontcam_timevect_lessthan = epocs.FrontCam_start_time(epocs.FrontCam_start_time <= time_vect_end); 
                FrontCamOnLocations = zeros(1,length(frontcam_timevect_lessthan)); %time vect is already downsampled
                for i = 1:length(frontcam_timevect_lessthan)
                    if epocs.FrontCam_start_time(i) <= time_vect_end 
                        tmp = find(time_vect >= epocs.FrontCam_start_time(i));  
                        FrontCamOnLocations(i) = tmp(1);  
                    end
                end

                epocs.FrontCamTTL = zeros(1,length(time_vect));
                for i=1:length(FrontCamOnLocations)
                    epocs.FrontCamTTL(FrontCamOnLocations(i)) = 5; 
                end
                clear frontcam_timevect_lessthan
            end
                
            %% Upsample the body position (Anymaze)
            if BehaviorAnymaze == 1
                if ~isempty(Body_parts)
                    for i = 1:length(Body_parts)
                        resample_Position.(AnimalID).(Body_parts{i}) = nan(1,size(BehavData.(AnimalID).Speed,2));   % inpaint_nans function allows to upsample
                    end
                    for i = 1:length(Body_parts)
                        resample_Position.(AnimalID).(Body_parts{i})(1,:) = inpaint_nans(BehavData.(AnimalID).(Body_parts{i})(1,:));   % inpaint_nans function allows to upsample
                    end
                end
            end
            
            %% Upsample the body position (Deeplabcut)
            if BehaviorDeeplabcut == 1
                if AnalyzeDeepLabCut == 1
                resample_Position = cell(2,1);
                for cam = 1:length(cameras)
                    if ~isempty(Position{cam})
                        tmpfld = fieldnames(Position{cam});
                        tmpfld2 = fieldnames(Position{cam}.(tmpfld{1}));
                        for i = 1:length(tmpfld)
                            for ii = 1:length(tmpfld2)
                %                 resample_Position{cam}.(tmpfld{i}).(tmpfld2{ii}) = nan(size(Position{cam}.(tmpfld{i}).(tmpfld2{ii})));
                                resample_Position{cam}.(tmpfld{i}).(tmpfld2{ii}) = inpaint_nans(Position{cam}.(tmpfld{i}).(tmpfld2{ii}));
                            end
                        end
                    end
                end
                end
            end
            
            %% Extract specific body part position: SPEED
            if BehaviorAnymaze == 1
                if ~isempty(Body_parts)
                    BodySpeed = resample_Position.(AnimalID).Speed;  
                %     [data_fp] = [data_fp;BodySpeed];

                    figure; clf
%                     if show_plot == 0
                           set(gcf,'visible','off')
%                     end
                    plot(time_vect,BehavData.(AnimalID).Speed,'.b') % raw data from Anymaze 
                    hold on
                    plot(time_vect,BodySpeed,'r') % upsampled data
                    yline(0.03)
                    yline(0.06)
                    min_x=-30; max_x=max(time_vect)+30; min_y=mean(BodySpeed)-0.2; max_y=mean(BodySpeed)+0.6; axis([min_x max_x min_y max_y]); 
                    xlabel('time(sec)'); ylabel('speed (m/s)'); title('Speed data'); 
                    L=legend('Speed upsampled','Speed raw'); L.Location = 'Best';

                    if  done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                        saveas(gcf,[PATH2SAVE,'figures\Speed upsampling, ',AnimalID,'.tif'])
                        saveas(gcf,[PATH2SAVE,'figures\Speed upsampling, ',AnimalID,'.fig'])
                    end
                end
            end
            % clear resample_Position field_names

            
              %% Extract specific body part position: SPEED
            if BehaviorDeeplabcut == 1
                if AnalyzeDeepLabCut == 1
                if ~isempty(Body_parts)
                    cam = 1;
                    BodySpeed = resample_Position{cam}.BodyCenter.speed;  

                    figure; clf
%                     if show_plot == 0
                           set(gcf,'visible','off')
%                     end
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
            end
            % clear resample_Position field_names
            
           
            
            %% Calculate angle
%             set(gcf,'visible','on')

            if BehaviorDeeplabcut == 1
                if cumanglecalc == 1
                cam = 1;
                neck_x = resample_Position{cam}.Nose.x; % could also do BodyCenter
                neck_y = resample_Position{cam}.Nose.x;  % one point every 10ms 
                tail_x = resample_Position{cam}.TailBase.x;
                tail_y = resample_Position{cam}.TailBase.y;
                [angle_deg_final] = rotations_function_forpaper(neck_x,neck_y,tail_x,tail_y,time_vect,fiber_side,dt_ds);
                for o = 1:size(epocs.PulseStart_time)
                    plot([epocs.PulseStart_time(o) epocs.PulseStart_time(o)],[-50 50],'m')
                    plot([epocs.PulseStop_time(o) epocs.PulseStop_time(o)],[-50 50],'c')
                end
                end
            end
            
           
            
            
            if calcdFF == 1
                %% LOW PASS FILTER OF FP DATA TO 1Hz
                if lowpassfiltering == 1
                    %filter option 2: butter
                    ftype = 'low';
                    n = 2; % 2nd order filter,
                    Wn = lowpassfreq/((sampling_rate_ds)/2); %lowpassfreq defined above
                    % 0.5 Hz = 2 sec ; 1 Hz = 1 sec ; 2 Hz = 0.5 sec ; 3 Hz = 0.33 sec
                    [a,b] = butter(n,Wn,ftype);
                    for channel=1:length(channel_names)
                        for colore=1:length(color_names)
                            streams.lowfilt.(channel_names{channel}).(color_names{colore}) = filtfilt(a,b,double(streams.rawdata.(channel_names{channel}).(color_names{colore})));
                            % plot
                            figure; clf; 
    %                         if show_plot == 0
                               set(gcf,'visible','off')
    %                         end
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
                    % % dFF synapse; % Note: Synapse uses an exponential smooth to estimate the F0 
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
                        % % dFF synapse; % Note: Synapse uses an exponential smooth to estimate the F0 
                        
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
    %                 if show_plot == 0
                       set(gcf,'visible','off')
    %                 end
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
                        streams.dFF_dtr.(dFF_names{d}) = detrend(streams.dFF.(dFF_names{d})); 
                        %plot
                        figure; clf
    %                     if show_plot == 0
                           set(gcf,'visible','off')
    %                     end
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
                        baseline_window = 10; %ADJUST
                        baseline_dt = 5; %Moving steps for the threshold 
                        percentile = 8; %8th percentile
                        % New time vector
                        temp_t = time_vect - time_vect(1); %creating new time vector set to starting at 0, incase time_vect weas not starting at 0
                        dummie = 1:length(temp_t); %generates index of same amount of samples 
                        dummie = dummie(temp_t >= baseline_window); %value when the time is beyond the window --> so that we can apply the moving window calculations
                        idx_baseline_window = dummie(1); %index of the timestamp as of which we can start the 60sec moving window
                        dummie = 1:length(temp_t);
                        dummie = dummie(temp_t >= baseline_dt);
                        idx_baseline_dt = dummie(1); %index of the timestamp as of which we start the 10sec moving window, i.e. baseline dt
                        clear dummie temp_t
                        % Index for the moving steps
                        ix1 = 1:idx_baseline_dt:length_data-idx_baseline_window; %
                           ix2 = idx_baseline_window:idx_baseline_dt:length_data; %index for the second bound of the window
                        ix = [ix1' ix2']; %first and last index you will check in each window (first colum onset, second colum offset) and the rows are the windows
                        if ix(end,2) < length_data
                            [ix] = [ix;[length_data-idx_baseline_window length_data]]; t
                        end
                        % Calculate 8th percentile across Pre and Post datasets
                        dFF_perc8.(dFF_names{d}) = ones(1,length_data)*nan; 
                        dFF_f0.(dFF_names{d}) = ones(1,length_data)*nan;
                        % calculate the entire 8th percentile (=baseline)
                        for o = 1:size(ix,1) %loop goes thru rows, equivalent of number of windows
                            idx = ix(o,1):ix(o,2); 
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
    %                     if show_plot == 0
                           set(gcf,'visible','off')
    %                     end
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
                    n = 2; % 2nd order filter
                    Wn = highpassfreq/((sampling_rate_ds)/2); %highpassfreq
                    [a,b] = butter(n,Wn,ftype);

                    for d=1:length(dFF_names)
                        streams.dFF_hp.(dFF_names{d}) = filtfilt(a,b,double(streams.dFF.(dFF_names{d})));
                        % plot
                        figure; clf
    %                     if show_plot == 0
                           set(gcf,'visible','off')
    %                     end
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
    %                 if show_plot == 0
                       set(gcf,'visible','off')
    %                 end
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

            end
            %% ONCE DONE WITH MODIFYNG dFF in streams, add SPEED (can also modify later)          
            if BehaviorAnymaze == 1
                for d=1:length(dFF_names)   
                    streams.speed.(dFF_names{d}) = BodySpeed;
                end
            end
            
            
               %% ONCE DONE WITH MODIFYNG dFF in streams, add SPEED (can also modify later) AND ROTATIONS         
            if BehaviorDeeplabcut == 1
                if AnalyzeDeepLabCut == 1
                for d=1:length(dFF_names)   
                    streams.speed.(dFF_names{d}) = BodySpeed;
                    streams.angle.(dFF_names{d}) = angle_deg_final';
                    streams.cumangle.(dFF_names{d}) = NaN(1,length(BodySpeed));
                end
                end
            end

            %% TRIGGER DATA BASED ON OPTOGENETIC STIMULATION (OR OTHER EVENT)
            % Setup the time vectors for each trial based on how big you want the window to be
            dummie = 1:length(time_vect); % indices 1 to the end of the time vector
            dummie = dummie(time_vect >= abs(TRANGE(1))); %select the indices of the time vector starting from the beginning of the time window
            idx_Init = dummie(1); % this was basically a trick to get the number of time vector indices needed to produce a duration of TRANGE(1) i.e. 15 sec
            dummie = 1:length(time_vect); 
            dummie = dummie(time_vect >= abs(TRANGE(2)));
            idx_End = dummie(1);
            n = idx_Init + idx_End; % Length of each trial in indices
            t_trials = time_vect(1:n) - time_vect(idx_Init); % get the time vector for 1 trial, of duration n, starting at initial index

            % Create 2-column array with the event to align to (start and stop)
            Events2Align2 = {'PulseStart_time','PulseStop_time'}            
            t_opto = [epocs.(Events2Align2{1}) epocs.(Events2Align2{2})]; % t_opto: a double column with all trials: column 1 (all rows) the start of the opto events and column 2 (all rows) the end
            if t_opto(end,2) == Inf %eg if trial was not finished
                t_opto(length(t_opto),:) = [];
            end

            % define the t_opto2 matrix that will contain the individual trials for each trial type (rows) - onset/offset (2 columns)
            % first look through the list of trials based on the power or frequency and group those that go together into a list of trial numbers 
            if ExperimentType == 'power'
                % % OPTION 1: POWER DEPENDENCE (CHRIMSON)
                if stim_duration == 10  % 10 seconds
                    pulse_number = 201
                elseif stim_duration == 3   % 3 seconds
                    pulse_number = 61
                else
                end
                StimAmpOne = epocs.StimAmp(1:pulse_number:length(epocs.StimAmp)-pulse_number-1); %isolate the first one of each series
                StimTimeOne = epocs.StimTime(1:pulse_number:length(epocs.StimTime)-pulse_number-1);
                PulseStart_type = zeros(size(epocs.PulseStart_time,1),1);
                for i=1:length(epocs.PulseStart_time)
                    for j=1:length(StimTimeOne)
                        if abs(epocs.PulseStart_time(i) - StimTimeOne(j)) < 0.1
                            PulseStart_type(i) = StimAmpOne(j)
                        else  
                        end
                    end
                end
                % trial allocation
                trial_select_820 = find(PulseStart_type>700)';
                trial_select_0 = find(PulseStart_type==0)';
                trial_select_160 = find(PulseStart_type>0 & PulseStart_type<300)';
                % min length
                minlength = min([length(trial_select_0),length(trial_select_160),length(trial_select_820)])
                % % % minlength = 12;
                % final trial allocation after correcting min length
                trial_select_820 = trial_select_820(1:minlength);
                trial_select_160 = trial_select_160(1:minlength);
                trial_select_0 = trial_select_0(1:minlength);            
                for d=1:length(dFF_names)
                    t_opto2.(dFF_names{d}).('p0uW') = t_opto(trial_select_0,:);  
                    t_opto2.(dFF_names{d}).('p500uW') = t_opto(trial_select_160,:);
                    t_opto2.(dFF_names{d}).('p2000uW') = t_opto(trial_select_820,:);
                end
                
            elseif ExperimentType == 'poLED'
                if contains(Virus,'465')
                    t_opto2.(dFF_names{d}).('p0uW') = t_opto([1,2,3,4],:); 
                    t_opto2.(dFF_names{d}).('p30uW') = t_opto([5,6,7,8],:);
                    t_opto2.(dFF_names{d}).('p350uW') = t_opto([21,22,23,24],:);
                end
                
            elseif ExperimentType == 'freq_'
            % %% OPTION 2: FREQUENCY DEPENDENCE normal
                PulseFreq = ones(size(epocs.PulseStart_time,1),1)*nan;
                idx_vect = ones(size(epocs.PulseStart_time,1),1)*nan;
                for i = 1:length(epocs.PulseStart_time)
                        temp = epocs.StimTime - epocs.PulseStart_time(i);
                        temp2=find(temp>0)
                        temp3 = min(temp2)
                        idx_vect(i) = temp3(1)
                end
                for i=2:length(idx_vect)
                    if idx_vect(i)- idx_vect(i-1) == 1
                        PulseFreq(i-1) = 0
                    elseif idx_vect(i)- idx_vect(i-1) == 101
                        PulseFreq(i-1) = 10
                    elseif idx_vect(i)- idx_vect(i-1) == 201
                        PulseFreq(i-1) = 20
                    elseif idx_vect(i)- idx_vect(i-1) == 401
                        PulseFreq(i-1) = 40
                    else
                        PulseFreq(i-1) = NaN;
                    end
                end
                if length(epocs.StimTime) - idx_vect(end) + 1 == 1
                    PulseFreq(length(idx_vect)) = 0
                elseif length(epocs.StimTime) - idx_vect(end) + 1 == 101
                    PulseFreq(length(idx_vect)) = 10
                elseif length(epocs.StimTime) - idx_vect(end) + 1 == 201
                    PulseFreq(length(idx_vect)) = 20
                elseif length(epocs.StimTime) - idx_vect(end) + 1 == 401
                    PulseFreq(length(idx_vect)) = 40
                else
                end
                % trial allocation
                trial_select_20 = find(PulseFreq==20)';
                trial_select_10 = find(PulseFreq==10)';
                trial_select_0 = find(PulseFreq==0)';
                minlength = min([length(trial_select_0),length(trial_select_10),length(trial_select_20)])
                  trial_select_20 = trial_select_20(1:minlength);
                trial_select_10 = trial_select_10(1:minlength);
                trial_select_0 = trial_select_0(1:minlength);       
                for d=1:length(dFF_names)
                    t_opto2.(dFF_names{d}).('p0Hz') = t_opto(trial_select_0,:); 
                    t_opto2.(dFF_names{d}).('p10Hz') = t_opto(trial_select_10,:);
                    t_opto2.(dFF_names{d}).('p20Hz') = t_opto(trial_select_20,:);
                end

            elseif ExperimentType == 'mobil' % 
                if stim_duration == 10
                    pulse_number = 201;
                elseif stim_duration == 3
                    pulse_number = 61;
                elseif stim_duration == 5
                    pulse_number = 101;
                else
                end
                % Identify Pulse start type
                StimAmpOne = epocs.StimAmp(1:pulse_number:length(epocs.StimAmp)-pulse_number-1); %isolate the first one of each series
                StimTimeOne = epocs.StimTime(1:pulse_number:length(epocs.StimTime)-pulse_number-1);
                PulseStart_type = zeros(size(epocs.PulseStart_time,1),1);
                for i=1:length(epocs.PulseStart_time)
                    for j=1:length(StimTimeOne)
                        if abs(epocs.PulseStart_time(i) - StimTimeOne(j)) < 0.1
                            PulseStart_type(i) = StimAmpOne(j);
                        else  
                        end
                    end
                end
                % Trial allocation
                trial_select_820 = find(PulseStart_type~=0)';    
                trial_select_0 = find(PulseStart_type==0)';
                minlength = min([length(trial_select_0),length(trial_select_820)])
                % minlength = 12;
                trial_select_820 = trial_select_820(1:minlength);
                trial_select_0 = trial_select_0(1:minlength);
                % Input trials into t_opto2
                t_opto2.(dFF_names{d}).('p0uW') = t_opto(trial_select_0,:); % t_opto2: group the ones to average from t_opto 
                if anymazepower == 2
                    t_opto2.(dFF_names{d}).('p2000uW') = t_opto(trial_select_820,:);
                elseif anymazepower == 0.2
                    t_opto2.(dFF_names{d}).('p200uW') = t_opto(trial_select_820,:);
                end
            
            end


            
            %% Generate the split dFF data (and other datatypes) trial matrixes based on trial allocation (t_opto2)
            for d=1:length(dFF_names)
                % Compute the aligned results of the trials, initialization to have them the number of rows that are in t_opto2 (=trials) and the length of the window: n
                for pow = 1:length(Opto_Powers)
                    for j = 1:length(datatype)
                        IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{j}) = ones(size(t_opto2.(dFF_names{d}).(Opto_Powers{pow}),1),n)*nan;
                        IndivStim_data.baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{j}) = ones(size(t_opto2.(dFF_names{d}).(Opto_Powers{pow}),1),n)*nan;
                    end
                end
                raw_or_corr = {'raw','baselinecorr'};
                
                % Compute the aligned results of the trials, assign the real data
                for pow = 1:length(Opto_Powers)
                    for i = 1:length(datatype)
                        for o = 1:size(t_opto2.(dFF_names{d}).(Opto_Powers{pow}),1) 
                            ix = find(abs(time_vect-t_opto2.(dFF_names{d}).(Opto_Powers{pow})(o,1)) == min(abs(time_vect-t_opto2.(dFF_names{d}).(Opto_Powers{pow})(o,1))));
                            if ix+idx_End <= length(time_vect) && ix-(idx_Init-1) >= 1 
                                tmp = ix - (idx_Init-1):ix + idx_End;
                                IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{i})(o,:) = streams.(datatype{i}).(dFF_names{d})(tmp);
                            elseif ix+idx_End > length(time_vect) 
                                tmp = ix - (idx_Init-1):length(time_vect); 
                                IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{i})(o,1:length(tmp)) = ...
                                        streams.(datatype{i}).(dFF_names{d})(tmp);
                            elseif ix-(idx_Init-1) < 1 
                                tmp = 1:(ix + idx_End); 
                                IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{i})(o,(n-length(tmp)+1):end) = ...
                                        streams.(datatype{i}).(dFF_names{d})(tmp);
                            end
                        end
                        % baseline correction 
                        IndivStim_data.baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{i})= ...
                        IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{i}) - nanmedian...
                        (IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{i})...
                        (:,t_trials >= BASELINE_WIN.OptoStim(1) & t_trials <= BASELINE_WIN.OptoStim(2)),2); %t_trials is the time vector for each trial
                     end
                end
            end

            
           
            
            %% CUMULATIVE ANGLE --> FROM STIM ONSET
            
            fromstimonset = 1; %if you want the cumulative angle calc to start only at stim onset
%             fromstimonset = 0; % 0 if you want the cumulative angle calc to start at the -15 sec or in the pre period
            
            if BehaviorDeeplabcut == 1
                if cumanglecalc == 1
                if fromstimonset == 0
                    for d=1:length(dFF_names)
                        for pow = 1:length(Opto_Powers)
                            for r = 1:length(raw_or_corr)
                                IndivStim_data.(raw_or_corr{r}).(dFF_names{d}).(Opto_Powers{pow}).cumangle = ...
                                    cumsum(IndivStim_data.(raw_or_corr{r}).(dFF_names{d}).(Opto_Powers{pow}).angle,2,'omitnan');                       
                                IndivStim_data.(raw_or_corr{r}).(dFF_names{d}).(Opto_Powers{pow}).cumangle = ...
                                    cumsum(IndivStim_data.(raw_or_corr{r}).(dFF_names{d}).(Opto_Powers{pow}).angle,2,'omitnan');
                            end
                        end
                    end
                
                    
                    
                elseif fromstimonset == 1
                    clear cumangle_full cumangle
                    i=3; % angle     
                    for d=1:length(dFF_names)
                        for pow = 1:length(Opto_Powers)
                            for r = 1:length(raw_or_corr)
                                for w=1:size(IndivStim_data.(raw_or_corr{r}).(dFF_names{1}).(Opto_Powers{pow}).angle,1)
                                    optostart = find(t_trials ==0);
                                    data2process = IndivStim_data.(raw_or_corr{r}).(dFF_names{1}).(Opto_Powers{pow}).(datatype{i})(w,optostart:size(IndivStim_data.(raw_or_corr{r}).(dFF_names{1}).(Opto_Powers{pow}).angle,2)); % 
                                    length_data = length(data2process);
                                    cumangle.(Opto_Powers{pow})(w,:) = ones(1,length_data)*nan;  
                                    cumangle.(Opto_Powers{pow})(w,:) = cumsum(data2process,2,'omitnan');                                               
                                end

                            cumangle_full.(raw_or_corr{r}).(Opto_Powers{pow}) = ones(size(IndivStim_data.raw.(dFF_names{1}).(Opto_Powers{pow}).angle,1),length(t_trials))*nan;  % first take vector with NaN, so that everything before stim onset becomes NaNs
                            cumangle_full.(raw_or_corr{r}).(Opto_Powers{pow})(:,optostart:end) = cumangle.(Opto_Powers{pow});
                            end
                            IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).cumangle = cumangle_full.raw.(Opto_Powers{pow});
                            IndivStim_data.baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).cumangle = cumangle_full.baselinecorr.(Opto_Powers{pow});
                        end
                    end
                   
                end
            end
             
            % baseline corr cum angle again
            for pow = 1:length(Opto_Powers)                        
                    if cumanglecalc == 1
                IndivStim_data.baselinecorr.(dFF_names{1}).(Opto_Powers{pow}).cumangle = ...
                    IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).cumangle - nanmedian...
                        (IndivStim_data.raw.(dFF_names{1}).(Opto_Powers{pow}).cumangle ...
                (:,t_trials >= BASELINE_WIN.OptoStim(1) & t_trials <= BASELINE_WIN.OptoStim(2)),2); %t_trials is the time vector for each trial
                
                IndivStim_data.baselinecorr.(dFF_names{1}).(Opto_Powers{pow}).cumangle = ...
                    IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).cumangle - nanmedian...
                        (IndivStim_data.raw.(dFF_names{1}).(Opto_Powers{pow}).cumangle ...
                (:,t_trials >= BASELINE_WIN.OptoStim(1) & t_trials <= BASELINE_WIN.OptoStim(2)),2); %t_trials is the time vector for each trial 
                end
            
            end
            
            
            
            end
   
            %%
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
            %% Plot the results per trial types separately
            wanttrialsepgraphs = 0;
            if wanttrialsepgraphs ==1
            trial_mode = {'raw','baseline_corrected'}; % raw or baseline_corrected.
            
            for d=1:length(dFF_names)
                % Individual trials
                avg_mode = 1; % 1, use mean and SEM; 2, use median and MAD
                for i = 1:length(datatype)
                    for r = 1:length(raw_or_corr)
                        figure;        
                        if show_plot == 0
                           set(gcf,'visible','off')
                        end
%                         scrsz = get(0,'ScreenSize'); %if you want screen size graphs
%                         set(gcf,'Position',scrsz);
                        for pow = 1:length(Opto_Powers)
                            subplot(length(Opto_Powers),1,pow); 
                            plot(t_trials,IndivStim_data.(raw_or_corr{r}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{i})); % plot individual trials
                            hold on
                            if avg_mode == 1
                                plot(t_trials,nanmean(IndivStim_data.(raw_or_corr{r}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{i}),1),'b',...
                                    'LineWidth',1.5) % plot the average on top of the individual trials
                            elseif avg_mode == 2
                                plot(t_trials,nanmedian(IndivStim_data.(raw_or_corr{r}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{i}),1),'b',...
                                    'LineWidth',1.5)
                            end
                            xline(0,'-k');
                            xline(stim_duration,'-m');
                            yline(0,'-.k');
                            if ~isempty(stim_duration2)
                                xline(stim_duration2,'-c');
                            end
                            xlabel('Time (s)')
                            ylabel(datatype{i});
                            if i == 4; % for speed
                                ylabel([datatype{i},' m/s']);
                            end
                            xlim([TRANGE(1) TRANGE(2)])
%                             ylim([round(min(streams.(datatype{i}).(dFF_names{d})))-5 round(max(streams.(datatype{i}).(dFF_names{d})))+5])
                            if i == 4; % for speed
                                ylim([round(min(streams.(datatype{i}).(dFF_names{d})))-0.1 round(max(streams.(datatype{i}).(dFF_names{d})))+0.5])
                            end
                            title([Opto_Powers{pow}(2:end)])
                        end
                        sgtitle([raw_or_corr{r},' ',datatype{i},' ',dFF_names{d}]); % for groups of subplots
                        if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                            saveas(gcf,[PATH2SAVE,'figures\INDIV TRIALS ',raw_or_corr{r},' ',datatype{i},' ',dFF_names{d},'.tif'])
                            saveas(gcf,[PATH2SAVE,'figures\INDIV TRIALS ',raw_or_corr{r},' ',datatype{i},' ',dFF_names{d},'.fig'])
                        end
                    end
                end

                % Trials, mean and SEM    
                color2plot = {'b','g','m','r','y','k'};
                % color2plot = {'r'};
                for i = 1:length(datatype)
                    for r = 1:length(raw_or_corr)
                        figure;         
                        if show_plot == 0
                           set(gcf,'visible','off')
                        end
%                         scrsz = get(0,'ScreenSize'); %if you want screen size graphs
%                         set(gcf,'Position',scrsz);
                        for pow = 1:length(Opto_Powers)
                            subplot(length(Opto_Powers),1,pow);
                            tmp_avg = nanmean(IndivStim_data.(raw_or_corr{r}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{i}),1); % plot the mean
                            tmp_error = nanstd(IndivStim_data.(raw_or_corr{r}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{i}),1,1)./...
                                sqrt(sum(~isnan(IndivStim_data.(raw_or_corr{r}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{i})(:,1)))); % plot the error
                            error_area(t_trials,tmp_avg,tmp_error,color2plot{pow},0.25)
                            xline(0,'-k');
                            xline(stim_duration,'-m');
                            if ~isempty(stim_duration2)
                                xline(stim_duration2,'-c');
                            end
                            yline(0,'-.k');
                            xlabel('Time (s)')
                            ylabel(datatype{i});
                            if i == 4; % for cum angle
                                ylabel([datatype{i},' (deg)']);
                            end
                            if i == 2; % for speed
                                ylabel([datatype{i},' cm/s']);
                            end
                            xlim([TRANGE(1) TRANGE(2)])
%                             ylim([round(min(streams.(datatype{i}).(dFF_names{d}))*1.4)-1 round(max(streams.(datatype{i}).(dFF_names{d}))*1.2)])
                            if i == 4; % for speed
                                ylim([round(min(streams.(datatype{i}).(dFF_names{d})))-0.1 round(max(streams.(datatype{i}).(dFF_names{d})))+0.4])
                            end
                            title([Opto_Powers{pow}(2:end)])
                        end
                        sgtitle([raw_or_corr{r},' ',datatype{i},' ',dFF_names{d}]); % for groups of subplots
                        if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                            saveas(gcf,[PATH2SAVE,'figures\AVE INDIV TRIALS RAW ',raw_or_corr{r},' ',datatype{i},' ',dFF_names{d},'.tif'])
                            saveas(gcf,[PATH2SAVE,'figures\AVE INDIV TRIALS RAW ',raw_or_corr{r},' ',datatype{i},' ',dFF_names{d},'.fig'])
                        end
                    end
                end
            end %for d=1:length(dFFnames)
            end
            %% Average plot
            
%             set(gcf,'visible','on')

            
            for d=1:length(dFF_names)
                %RAW
                color2plot = {'b','g','m','r','c','k'};
                for i=1:length(datatype)
                    for r=1:length(raw_or_corr)
                        figure; 
                        if show_plot == 0
                            set(gcf,'visible','off')
                        end
                        hold on
                        for pow = 1:length(Opto_Powers)
                            tmp_avg = nanmean(IndivStim_data.(raw_or_corr{r}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{i}),1);
                            tmp_error = nanstd(IndivStim_data.(raw_or_corr{r}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{i}),1,1)./...
                                sqrt(sum(~isnan(IndivStim_data.(raw_or_corr{r}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{i})(:,1))));
                            error_area(t_trials,tmp_avg,tmp_error,color2plot{pow},0.25); % error_area(X, Y, barlength, color, alpha, varargin) to plot SEM
%      
                        end
                        xline(0,'-k'); %other option: plot([0 0],limits2plot.dLight,'k')
                        xline(stim_duration,'-m');
                        yline(0,'-.k');
%                         yline(max(max_avg),'-.r');
                        xlabel('Time (s)');
                        ylabel(datatype{i});
                        if i == 4; % for cum angle
                                ylabel([datatype{i},' (deg)']);
                        end
                        if i == 2; % for speed
                            ylabel([datatype{i},' cm/s']);
                        end
                        xlim([TRANGE(1) TRANGE(2)]);
%                         ylim(limits2plot.(raw_or_corr{r}).(dFF_names{d}).(datatype{i})); %see parameters setting at start                     
%                         PeakAve = num2str(mean(max_avg));
%                         PeakMax = num2str(max(max_avg));
                        sgtitle([raw_or_corr{r},' ',datatype{i},' ',dFF_names{d}],'Interpreter','none')
%                         sgtitle([raw_or_corr{r},' ',datatype{i},' ',dFF_names{d},', AvePeak: ',PeakAve,', MaxPeak: ',PeakMax],'Interpreter','none')

                        % legend
                        if TrialType_Number == 2
                        % 2 powers
                        legend([Opto_Powers{1}(2:end),'error'],[Opto_Powers{1}(2:end)],...
                                [Opto_Powers{2}(2:end),'error'],[Opto_Powers{2}(2:end)],...
                                'Opto On','Opto off','Location','northwest','NumColumns',1)      
                        elseif TrialType_Number == 3
                        % 3 powers
                        legend([Opto_Powers{1}(2:end),'error'],[Opto_Powers{1}(2:end)],...
                            [Opto_Powers{2}(2:end),'error'],[Opto_Powers{2}(2:end)],...
                            [Opto_Powers{3}(2:end),'error'],[Opto_Powers{3}(2:end)],...
                            'Opto On','Opto off','Location','northwest','NumColumns',1)
                        elseif TrialType_Number == 4
                        % 4 powers
                        legend([Opto_Powers{1}(2:end),'error'],[Opto_Powers{1}(2:end)],...
                            [Opto_Powers{2}(2:end),'error'],[Opto_Powers{2}(2:end)],...
                            [Opto_Powers{3}(2:end),'error'],[Opto_Powers{3}(2:end)],...
                            [Opto_Powers{4}(2:end),'error'],[Opto_Powers{4}(2:end)],...
                            'Opto On','Opto off','Location','northwest','NumColumns',1)
                        elseif TrialType_Number == 5
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
                            saveas(gcf,[PATH2SAVE,'figures\AVE all powers optostim ',raw_or_corr{r},' ',datatype{i},' ',dFF_names{d},'.tif'])
                            saveas(gcf,[PATH2SAVE,'figures\AVE all powers optostim ',raw_or_corr{r},' ',datatype{i},' ',dFF_names{d},'.fig'])
                        end
                    end
                end
            end %for d=1:length(dFFnames)
            
        

            
        end % end of if loop "if done" --> maybe need to put this a bit lower
 
      
        
        
        %% SAVE INDIVDIDUAL SESSION    --> ADD TO SAVE ADDITIONAL THINGS LIKE CORRELATIONS LATER
        
        if BehaviorAnymaze == 0
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
            save([PATH2SAVE,'IndividualData.mat'],'Params','IndivStim_data','resample_Position','epocs','t_trials','datatype','dFF_names','Opto_Powers','Events2Align2','BASELINE_WIN',...
                'stim_duration','TRANGE','dt_ds','sampling_rate_ds','length_data','time_vect','streams');
            % NOT SURE WE NEED streams and time_vect
        end    


        %% Ending all the loops

    end %for s=1:length sessions
end % for all mice: i=1:numfiles, see above

    
