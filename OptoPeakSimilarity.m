% OptoPeakSimilarity.m
% Matlab code for new analysis for reanalysis of ondividual opto trials on PFC dLight Opto data
% March 2024

% need function 'error_area'

%% INITIALIZATIONS
close all; clear variables; clc;
set(0,'defaultfigurecolor',[1 1 1])

%% PARAMETERS 
show_plot = 1; % If 0, plots are not displayed
save_plot = 1; % If 0, plots are not saved
reanalysis = 0; % If 1, the code runs for sessions already analyzed. If 0, session excluded if analysed. This is based on the existence of the matlab space "IndividualData.mat"
overwrite = 0; % If reanalyzing data, 1 if you want to save them and overwrite results and 0 if not          
pooledtype = {'raw','baselinecorr'};

%% Define the path where the data is
PATH2DATA_0 = uigetdir('select folder'); %select the folder above the matlab space: Z:\Wolfrum-Group\Marie\UZH\DATA_PFC_dLight_Opto_DoseResponse\PFCdLightOPTO\_ALLDATA\DETQ10
% TrialDay = {'VEH','DETQ'}; 
PATH2SAVEPOOL = PATH2DATA_0;
mkdir([PATH2SAVEPOOL,'\pooled figures\']);
% Opto_Powers = {'SAL','COC'}; %  % imported via matlab .mat


%% Load data
load([PATH2DATA_0,'\PooledAllMice_DETQ10.mat']);


%% Plot
figure; 

% ID 0154
for i=1:14 % select 14 out of 18 trials
    plot(PooledINDIV_VEHnormalized.DETQ.baselinecorr.PFC.DRUG.dFFwithin.ID_0154.SessionNum1(i,:)); hold on; 
    title(['ID 0154'])
end


%% Downsample
trialnumb = {'t1','t2','t3','t4','t5','t6','t7','t8','t9','t10','t11','t12','t13','t14'}

downsample_freq = 0.4;
downsample_interval = 1/downsample_freq;
%1fps --> 120 points, 1 per second, i.e. 60 in the opto phase. %0.5fps --> 60 points, 1 per 2 second, i.e. 30 in the opto phase.

for ttt = 1:length(trialnumb)
    for nummice = 1:length(mice_list_virus.DETQ)
        data2plot=PooledINDIV_VEHnormalized.DETQ.baselinecorr.PFC.DRUG.dFFwithin.(PooledAnimalID_only.DETQ{nummice}).SessionNum1(ttt,:);
        datattt = timetable(data2plot','SampleRate',sampling_rate_ds);
        datatttAvg = retime(datattt,'regular','mean','SampleRate',downsample_freq); 
        PooledINDIV_VEHnormalized_ds.DETQ.baselinecorr.PFC.DRUG.dFFwithin.(PooledAnimalID_only.DETQ{nummice}).SessionNum1(ttt,:) = datatttAvg.Var1;
    end
                
    t_trials_ds = [seconds(datatttAvg.Time)] + t_trials(1); %also add the first timepoint from the previous timevector allowing to shift in time                
    t_trials_ds_minutes = t_trials_ds/60;
                

end


%% Calculate Average and SEM
    color2plot = {'b','g','r','k','c','m','y','b','g','r','k','c','m','y','b','g','r','k','c','m','y'};

% initialize. Use dFF_within and VEH_normalized
for ttt = 1:length(trialnumb)
    PooledINDIV_VEHnormalized_pool.(trialnumb{ttt}) = ones(length(mice_list_virus.DETQ),size(PooledINDIV_VEHnormalized_ds.DETQ.baselinecorr.PFC.DRUG.dFFwithin.ID_0154.SessionNum1(1,:),2))*nan;
end

% create a pool
for ttt = 1:length(trialnumb)
    for nummice = 1:length(mice_list_virus.DETQ)
        PooledINDIV_VEHnormalized_pool.(trialnumb{ttt})(nummice,:) = PooledINDIV_VEHnormalized_ds.DETQ.baselinecorr.PFC.DRUG.dFFwithin.(PooledAnimalID_only.DETQ{nummice}).SessionNum1(ttt,:);  
    end
end

% calculate AVE and SEM
for ttt = 1:length(trialnumb)
    Merged_INDIV_AVE(ttt,:) = nanmean(PooledINDIV_VEHnormalized_pool.(trialnumb{ttt}),1);  
    Merged_INDIV_SEM(ttt,:) = nanstd(PooledINDIV_VEHnormalized_pool.(trialnumb{ttt}),1,1)./sqrt(size(PooledINDIV_VEHnormalized_pool.(trialnumb{ttt}),1));  
end

% plot with SEM
figure;
for ttt = 1:length(trialnumb)
    error_area(t_trials_ds,Merged_INDIV_AVE(ttt,:),Merged_INDIV_SEM(ttt,:),color2plot{ttt},0.25); hold on; 
end

% plot without SEM
figure;
for ttt = 1:length(trialnumb)
    plot(t_trials_ds,Merged_INDIV_AVE(ttt,:),color2plot{ttt}); hold on; 
end
xline(0); xline(2)

%% Create a table for direct copy-paste into GraphPad

Data4GraphPad = ones(size(PooledINDIV_VEHnormalized_pool.(trialnumb{ttt}),2),5*length(trialnumb))*nan;

for ttt = 1:length(trialnumb)
    Data4GraphPad(:,ttt*5-5+1:ttt*5-5+1+4) = PooledINDIV_VEHnormalized_pool.(trialnumb{ttt})' 
end

 










