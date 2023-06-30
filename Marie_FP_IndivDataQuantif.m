% Marie_FP_IndivData_quantification
% Marie Labouesse, marie.labouesse@gmail.com - Feb 2021

% INDIVIDAL TRIALS AND MEAN, PEAK/DIP QUANTIFICATION
% load matlab spaces generated in: Marie_FP_IndivData_extraction and PooledData
% uses individual trials aligned to specific event to calculate basic quantifications to characterize the peak or dip: 
% quantifications are done on individual trials as well as on the mean of all trials for 1 mouse (in 1 session): AUC, min/max, latency, degree of inhibition (for example)
% generates and saves the data into matlab spaces and tables that can be copy-pasted into GraphPad for stats

% % INDIVIDAL TRIALS AND MEAN, CORRELATION DFF/BEHAVIOR QUANTIFICATION 
% calculates individual correlation between trial dFF quantifications and behavior for each trial
        % for now we want to correlate average speed "in the 5sec before stim winsow" with average Zscore dFF
% generates and saves the data into matlab spaces, graphs tables that can be copy-pasted into GraphPad for stats
% Also calculates the mean of correlation between trial dFF quantifications and behavior for each trial:
        % for now we want to correlate average speed "in the 5sec before stim winsow" with average Zscore dFF

% CONTROL: shuffle DATA
% calculate scrambled data (shuffled timeseries) to compare to this

%% INITIALIZATIONS
close all; clear variables; clc;
set(0,'defaultfigurecolor',[1 1 1])


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZATION

% LOAD MATLAB SPACES
PATH2DATA = uigetdir('select folder'); %path to folder above chrimson and mcherry
PATH2SAVE = PATH2DATA;
load([PATH2DATA,'\PooledAllMice.mat']);

% ANIMAL IDs
%LED465/595
% AnimalID_LED595=fieldnames(PooledAnimalID.LED595)'
% AnimalID_LED465=fieldnames(PooledAnimalID.LED465)'
%Chrimson/mCherry
AnimalID_Chrimson=fieldnames(PooledAnimalID.Chrimson)';
% AnimalID_mCherry=fieldnames(PooledAnimalID.mCherry)';

for v=1:length(virus)
    mice_list_virus.(virus{v}) = fieldnames(PooledAnimalID.(virus{v}))
end

% SESSIONS ..
SessionIDs = {'SessionNum1'}


   

%% Compute the AUC and the Peaks
% Determine the indexes where to calculate the AUC and other variables in the stimulation period
temp_t = time_vect - time_vect(1); %values in timevect shifted by one value, not sure why since time_vect(0) == . We dont use the t_trials cos it doesnt start at 0
dummie = 1:length(temp_t); % indexes of time vector
dummie = dummie(temp_t >= stim_duration); % so look for the indexes when the time vector is superior or equal to the stim duration
%                 dummie = dummie(temp_t >= stim_dur_per_block(1));
idx_AUC2 = dummie(1); %find the first index when the time vector hit the exact value of the end of the stim duration; this is the second AUC index, or the number of indixes you need to go thru to have the stim duration
idx_AUC1 = find(t_trials == 0); % this is the index in t_trials where the stimulation starts
idx_AUC = idx_AUC1:(idx_AUC1+(idx_AUC2-1)); % we will be calculating the AUC from the beginning of the stimulation to the index just before the first index at the end of the stim
clear dummie idx_AUC1 idx_AUC2
        
% Determine the indexes where to calculate the thresh period (5 seconds before stim onset)
thresh_duration = -5;
dummie = 1:length(t_trials); % indexes of time vector
dummie = dummie(t_trials >= thresh_duration);
idx_AUCthresh1 = dummie(1);
idx_AUCthresh2 = find(t_trials == 0); % this is the index in t_trials where the stimulation starts
idx_AUC_thresh = idx_AUCthresh1:idx_AUCthresh2;

% Determine the indexes where to calculate the pre period (10 seconds before the thresh onset, ie 15 seconds before the stim onset)
idx_AUCpre1 = 1; % ie first point of recording. 
pre_end = -5;
dummie = 1:length(t_trials); % indexes of time vector
dummie = dummie(t_trials >= pre_end);
idx_AUCthresh1 = dummie(1);
idx_AUCpre2 = idx_AUCthresh1-1; % this is the index in t_trials where the stimulation starts
idx_AUC_pre = idx_AUCpre1:idx_AUCpre2;



%% Initialize the "Measurements" structure where you will be putting the AUC, Peaks and specific values. 

for v=1:length(virus) %Chrimson or mCherry
    for p = 1:length(pooledtype)   % raw or baselinecorr
        for d=1:length(dFF_names)   % GPe and SNr (for example)
            for pow = 1:length(Opto_Powers) % 0 and 2mW (for example)
                for k = 1:size(datatype,2) % ZScoredFF and dFF and dFFsy and speed (for example) ---> NB might to analyze speed differently
                    for nummice=1:length(mice_list_virus.(virus{v}))    % all mice
                        for s = 1:length(SessionIDs) % all sessions for this mouse, eg SessionNum1 
                            t_opto_i = PooledINDIV.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(mice_list_virus.(virus{v}){nummice}).(SessionIDs{s});

                            Measurements.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(mice_list_virus.(virus{v}){nummice}).(SessionIDs{s}).AUC = ones(size(t_opto_i,1),1)*nan; %same number of rows as trials
                            Measurements.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(mice_list_virus.(virus{v}){nummice}).(SessionIDs{s}).DegInhibition = ones(size(t_opto_i,1),3)*nan; %value at start, end and degree of inhibition
                            Measurements.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(mice_list_virus.(virus{v}){nummice}).(SessionIDs{s}).LocalMinima = ones(size(t_opto_i,1),2)*nan; %value at minima, time at minima
                            Measurements.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(mice_list_virus.(virus{v}){nummice}).(SessionIDs{s}).LocalMaxima = ones(size(t_opto_i,1),2)*nan;
                            Measurements.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(mice_list_virus.(virus{v}){nummice}).(SessionIDs{s}).DegChange_Pre = ones(size(t_opto_i,1),3)*nan; %value at start, end and degree of inhibition
                            Measurements.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(mice_list_virus.(virus{v}){nummice}).(SessionIDs{s}).Average_Pre = ones(size(t_opto_i,1),1)*nan;
                            Measurements.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(mice_list_virus.(virus{v}){nummice}).(SessionIDs{s}).Average_Thresh = ones(size(t_opto_i,1),1)*nan;                            
                            Measurements.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(mice_list_virus.(virus{v}){nummice}).(SessionIDs{s}).Average_Post = ones(size(t_opto_i,1),1)*nan; 
                        
                        end
                    end
                end
            end
        end
    end
end

%% Determine the AUC for entire section (anything below 0 will count for negative AUC) 
% + the value at start of stim 
% + the value at end of stim 
% + the minimal value in the section
% + the local minima in section; by looking at when first derivative = 0
% + first inflection point after stim, ie when first derivative becomes negative = when inhibition starts.

% Divide the trace in each sector above or below 0. Get the AUC and peak for each section 
for v = 1:length(virus) %Chrimson or mCherry
    for p = 1:length(pooledtype)   % raw or baselinecorr
        for d = 1:length(dFF_names)   % GPe and SNr (for example)
            for pow = 1:length(Opto_Powers) % 0 and 2mW (for example)
                for k = 1:size(datatype,2) % ZScoredFF and dFF and dFFsy and speed (for example) ---> NB might to analyze speed differently
                    for nummice = 1:length(mice_list_virus.(virus{v}))    % all mice
                        for s = 1:length(SessionIDs)
                            t_opto_i = PooledINDIV.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(mice_list_virus.(virus{v}){nummice}).(SessionIDs{s});
                            fulldata = fullstreams.(virus{v}).(mice_list_virus.(virus{v}){nummice}).(SessionIDs{s}).(datatype{k}).(dFF_names{d});
                            firstderivative = gradient(fulldata);
                            for w = 1:size(t_opto_i,1)   % = how many rows
                                %AUC, includes positive and negative values
                                dff_all = t_opto_i(w,:); % put the trace of this row in a new variable cos easier to handle
                                dff_stim = t_opto_i(w,idx_AUC);   %piece of dFF in the stim epoch we are interested in for this row
                                tmp_AUC = trapz(dff_stim);
                                average_post = mean(dff_stim);
                                Measurements.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(mice_list_virus.(virus{v}){nummice}).(SessionIDs{s}).AUC(w) = tmp_AUC; %same number of rows as trials
                                %Value at start and end and amplitude of inhibition
                                dFFonset = t_opto_i(w,idx_AUC(1));
                                dFFoffset = t_opto_i(w,idx_AUC(end));
                                DegInh = dFFoffset - dFFonset; 
                                Measurements.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(mice_list_virus.(virus{v}){nummice}).(SessionIDs{s}).DegInhibition(w,:) = [dFFonset dFFoffset DegInh]; %value at start, end and degree of inhibition
                                % Same thing but for the PRE period
                                dFFonset_pre = t_opto_i(w,idx_AUC_pre(1));
                                dFFoffset_pre = t_opto_i(w,idx_AUC_pre(end));
                                DegInh_pre = dFFoffset_pre - dFFonset_pre; 
                                Measurements.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(mice_list_virus.(virus{v}){nummice}).(SessionIDs{s}).DegChange_Pre(w,:) = [dFFonset_pre dFFoffset_pre DegInh_pre]; %value at start, end and degree of inhibition

                                % Average Pre and Thresh periods
                                dff_stim_pre = t_opto_i(w,idx_AUC_pre);   %piece of dFF in the stim epoch we are interested in for this row
                                average_pre = mean(dff_stim_pre);
                                dff_stim_thresh = t_opto_i(w,idx_AUC_thresh);   %piece of dFF in the stim epoch we are interested in for this row
                                average_thresh = mean(dff_stim_thresh);
                                Measurements.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(mice_list_virus.(virus{v}){nummice}).(SessionIDs{s}).Average_Pre(w,:) = [average_pre]; %value at start, end and degree of inhibition
                                Measurements.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(mice_list_virus.(virus{v}){nummice}).(SessionIDs{s}).Average_Post(w,:) = [average_post]; %value at start, end and degree of inhibition
                                Measurements.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(mice_list_virus.(virus{v}){nummice}).(SessionIDs{s}).Average_Thresh(w,:) = [average_thresh]; %value at start, end and degree of inhibition
                                
                                
                                %Calculate first and second derivative
                                %first
                                dy = diff(dff_all); % 
                                dx = diff(t_trials); % 
                                dff_first = dy./dx; % 
                                t_first = 0.5*(t_trials(1:end-1)+t_trials(2:end)); 
                                %second
                                dyfirst = diff(dff_first); % 
                                dxfirst = diff(t_first); % 
                                dff_second = dyfirst./dxfirst; % 
                                t_second = 0.5*(t_first(1:end-1)+t_first(2:end)); 
%                                 plot(t_second,dff_second);
                                %Local Minima and Maxima in stim epoch and Latency to reach Minima and Maxima
                                      A = dff_first > 0; 
                                B = find(A == 1); 
                                C = find(diff(B) > 1);  
                                     positive_trans = ones(length(C),2)*nan;
                                negative_trans = ones(length(C),2)*nan;
                                for o = 1:length(C)
                                    positive_trans(o,:) = [B(C(o)) B(C(o))+1]; 
                                    negative_trans(o,:) = [B(C(o)+1)-1 B(C(o)+1)]; 
                                end
                                 %local minima
                                 for o = 1:length(C)
                                    localmaxima = max([dff_all(positive_trans(o,1):positive_trans(o,2)+1)]); 
                                    ix = find(dff_all(positive_trans(o,1):positive_trans(o,2)+1)==localmaxima);
                                    vectx = positive_trans(o,1):positive_trans(o,2)+1;
                                    latencytomaxima = t_trials(vectx(ix(1))); % 
                                    localminima = max([dff_all(negative_trans(o,1):negative_trans(o,2)+1)]);
                                    x = find(dff_all(negative_trans(o,1):negative_trans(o,2)+1)==localmaxima);
                                    vectx = negative_trans(o,1):negative_trans(o,2)+1;
                                    latencytominima = t_trials(vectx(ix(1))); % 
                                 end
                                 Measurements.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(mice_list_virus.(virus{v}){nummice}).(SessionIDs{s}).LocalMinima(w,:) = [localminima latencytominima];
                                 Measurements.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(mice_list_virus.(virus{v}){nummice}).(SessionIDs{s}).LocalMaxima(w,:) = [localmaxima latencytomaxima];                                
 
                                 
                            
                            end
                        end
                    end
                end
            end
        end
    end
end

                                
%% Generate the AVE Measurements data (average for each mouse, for all its trials) + easy "copy-paste" table for plotting into graph pad % then heatmaps (other script)
for v = 1:length(virus) %Chrimson or mCherry
    for p = 1:length(pooledtype)   % raw or baselinecorr
        for d = 1:length(dFF_names)   % GPe and SNr (for example)
            for pow = 1:length(Opto_Powers) % 0 and 2mW (for example)
                for k = 1:size(datatype,2) % ZScoredFF and dFF and dFFsy and speed (for example) ---> NB might to analyze speed differently
                    for nummice = 1:length(mice_list_virus.(virus{v}))    % all mice
                        for s = 1:length(SessionIDs)
                            if s == 1
                                Measurements_AVE.(pooledtype{p}).(dFF_names{d}).(datatype{k}).(virus{v}).AUC(nummice,pow) = ...
                                    nanmean(Measurements.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(mice_list_virus.(virus{v}){nummice}).(SessionIDs{s}).AUC,1)
                                Measurements_AVE.(pooledtype{p}).(dFF_names{d}).(datatype{k}).(virus{v}).DegInhibition(nummice,pow) = ...
                                    nanmean(Measurements.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(mice_list_virus.(virus{v}){nummice}).(SessionIDs{s}).DegInhibition(:,3),1)
                                Measurements_AVE.(pooledtype{p}).(dFF_names{d}).(datatype{k}).(virus{v}).DegChange_Pre(nummice,pow) = ...
                                    nanmean(Measurements.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(mice_list_virus.(virus{v}){nummice}).(SessionIDs{s}).DegChange_Pre(:,3),1) 
                                Measurements_AVE.(pooledtype{p}).(dFF_names{d}).(datatype{k}).(virus{v}).LocalMinima(nummice,pow) = ...
                                    nanmean(Measurements.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(mice_list_virus.(virus{v}){nummice}).(SessionIDs{s}).LocalMinima(:,1),1)
                                Measurements_AVE.(pooledtype{p}).(dFF_names{d}).(datatype{k}).(virus{v}).LocalMaxima(nummice,pow) = ...
                                    nanmean(Measurements.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(mice_list_virus.(virus{v}){nummice}).(SessionIDs{s}).LocalMaxima(:,1),1)
                                Measurements_AVE.(pooledtype{p}).(dFF_names{d}).(datatype{k}).(virus{v}).LatencyLocalMinima(nummice,pow) = ...
                                    nanmean(Measurements.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(mice_list_virus.(virus{v}){nummice}).(SessionIDs{s}).LocalMinima(:,2),1)
                                Measurements_AVE.(pooledtype{p}).(dFF_names{d}).(datatype{k}).(virus{v}).LatencyLocalMaxima(nummice,pow) = ...
                                    nanmean(Measurements.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(mice_list_virus.(virus{v}){nummice}).(SessionIDs{s}).LocalMaxima(:,2),1)                              
                                Measurements_AVE.(pooledtype{p}).(dFF_names{d}).(datatype{k}).(virus{v}).Average_Pre(nummice,pow) = ...
                                    nanmean(Measurements.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(mice_list_virus.(virus{v}){nummice}).(SessionIDs{s}).Average_Pre(:,1),1)
                                Measurements_AVE.(pooledtype{p}).(dFF_names{d}).(datatype{k}).(virus{v}).Average_Thresh(nummice,pow) = ...
                                    nanmean(Measurements.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(mice_list_virus.(virus{v}){nummice}).(SessionIDs{s}).Average_Thresh(:,1),1)
                                Measurements_AVE.(pooledtype{p}).(dFF_names{d}).(datatype{k}).(virus{v}).Average_Post(nummice,pow) = ...
                                    nanmean(Measurements.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(mice_list_virus.(virus{v}){nummice}).(SessionIDs{s}).Average_Post(:,1),1)  
                            elseif length(sessions) > 1
                            end
                        end
                    end
                end
            end
        end
    end
end



%% Save data
save([PATH2SAVE,'Cohort_analysis.mat'],'Measurements','Measurements_AVE','virus','PooledAnimalID','dFF_names','datatype','Opto_Powers',...
    't_trials','stim_duration','TRANGE','BASELINE_WIN','Events2Align2','dt_ds','sampling_rate_ds','length_data','fullstreams','time_vect','pooledtype');
    
    