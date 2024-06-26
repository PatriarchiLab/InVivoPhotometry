% ChemoSLA_PooledGraphs_2023
% Marie Labouesse, marie.labouesse@gmail.com - Sep 2023
% Cohort: NAc dLight1.3b weekly DETQ or food rest

% 1- POOLED DATA
% load matlab spaces generated in: WeeklyDETQ_SLA_DataExtraction
% puts individual trials into one big matrix (all mice together): PooledINDIV
% calculates the mean of individual trial trajectories aligned to specific event for each mouse: PooledAVE
% saves a workspace that can then be used for quantifications

% 2- POOLED GRAPHS
% generates average graphs curves (this used to be in the main script)

% NEED FUNCTIONS
%error_area_only_rectangle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
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
PATH2DATA_0 = uigetdir('select folder'); %select the overall folder, where the different groups are saved (as folders)- we will import all
virus = {'wk1','wk2','wk3'}; % weekly
% virus = {'FED','FOODRESTR','FED2'}; 


for v=1:length(virus)
    PATH2DATA.(virus{v}) = [PATH2DATA_0,'\',virus{v}];
    PATH2SAVEFOLDER.(virus{v}) = [PATH2DATA_0,'\',virus{v}];
    mice_list_virus.(virus{v}) = dir(PATH2DATA.(virus{v})); %all things in this folder
    mkdir([PATH2SAVEFOLDER.(virus{v}),'\pooled figures\']);
    mkdir([PATH2SAVEFOLDER.(virus{v}),'\pooled data\']);
    mkdir([PATH2DATA_0,'\pooled figures\']);
    mkdir([PATH2DATA_0,'\pooled data\']);
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

end

              

%% Import individual workspaces and load the data into a pooled array
for v=1:length(virus)
    %AnimalIDs
    for nummice=1:length(mice_list_virus.(virus{v}))
        AnimalIDs(nummice) = {['ID',mice_list_virus.(virus{v})(nummice).name(1:4)]};
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
        
        %% Loop for all the sessions for the mouse: Sessions are experimental days you want to average together (replicates of each other)- we only have 1 session here
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
            done = exist([PATH2SAVE,'PooledAllMice.mat'],'file'); % .................. 
            if done == 0 || reanalysis == 1 %if not done or if wish to reanalyze
        
                % load the mouse/session matlab space
                load([PATH2SESSION,'\IndividualData.mat']);
       
                % Initialization of pooled structure with individual trials   %................  ONLY ONE SESSION PER MOUSE.

                if nummice == 1 && s == 1
                    for d=1:length(dFF_names)
                        for pow = 1:length(Opto_Powers)
                            for k = 1:size(datatype,2)
                                % ANIMAL ID .................... add session ID later if needed
                                PooledAnimalID.(virus{v}) = {};
                                % STREAMS                
                                fullstreams.(virus{v}).(AnimalIDs{nummice}) = {};
                                % STREAMS ALIGNED                
                                fullstreams_aligned.(virus{v}).(Opto_Powers{pow}).(AnimalIDs{nummice}) = {};
                                % TIME VECT ALIGNED
                                timevect_aligned.(virus{v}).(Opto_Powers{pow}).(AnimalIDs{nummice}) = {}; 
                            end
                        end
                    end
                end
                
                for d=1:length(dFF_names)
                        for pow = 1:length(Opto_Powers)
                            for k = 1:size(datatype,2)
                                % Add data one animal at a time
                                % ANIMAL ID and SESSION ID
                                PooledAnimalID.(virus{v}) = AnimalIDs;

                                %STORE STREAMS
                                fullstreams.(virus{v}).(AnimalIDs{nummice}) = streams.dFF;

                                %STORE STREAMS ALIGNED TO OPTO ONSET #1
                                fullstreams_aligned.(virus{v}).(Opto_Powers{pow}).(AnimalIDs{nummice}) = streams_aligned.(datatype{k}).(dFF_names{d}).(Opto_Powers{pow});
                                timevect_aligned.(virus{v}).(Opto_Powers{pow}).(AnimalIDs{nummice}) = time_vect_aligned.(Opto_Powers{pow}); % DETQ/VEH or DRUG levels
                            end
                        end
                end
                
            end
        end
    end
end

%% Adjust all to VEH 0 as baseline
for v=1:length(virus)    
    for nummice=1:length(mice_list_virus.(virus{v}))
        for d=1 %:length(dFF_names)
            for k = 1 %:size(datatype,2)
                for pow = 1:length(Opto_Powers)
                    fullstreams_aligned_norm.(virus{v}).(Opto_Powers{pow}).(AnimalIDs{nummice}) = ...
                        fullstreams_aligned.(virus{v}).(Opto_Powers{pow}).(AnimalIDs{nummice}) ...
                        - mean(fullstreams_aligned.(virus{v}).VEH.(AnimalIDs{nummice}))                    
                end
            end
        end
    end
end

%if happy
fullstreams_aligned = fullstreams_aligned_norm;

%% Streams aligned pooled
   
    % create empty pooled array
    nummice = 1; k=2; d=1; pow = 1; s=1; 
    for v=1:length(virus)    
        for d=1:length(dFF_names)
            for pow = 1:length(Opto_Powers)
                for k = 1:size(datatype,2)
                    length_fullstreams_aligned.(Opto_Powers{pow}) = ...
                    size(find(~isnan(fullstreams_aligned.(virus{v}).(Opto_Powers{pow}).(AnimalIDs{nummice}))),2)       

                    fullstreams_aligned_pooled.(virus{v}).(Opto_Powers{pow}) = ...
                        ones(length(mice_list_virus.(virus{v})),length_fullstreams_aligned.(Opto_Powers{pow}))*nan;
                end
            end
        end

        %STORE STREAMS ALIGNED TO EVENT ONSET #1 with the average of each mouse in each row (eg 3 mice --> 3 rows)
        for nummice=1:length(mice_list_virus.(virus{v}))
            for d=1:length(dFF_names)
                for pow = 1:length(Opto_Powers)
                    for k = 1:size(datatype,2)
                        fullstreams_aligned_pooled.(virus{v}).(Opto_Powers{pow})(nummice,1:length_fullstreams_aligned.(Opto_Powers{pow})) = ...
                            fullstreams_aligned.(virus{v}).(Opto_Powers{pow}).(AnimalIDs{nummice})(1:length_fullstreams_aligned.(Opto_Powers{pow}));
                          
                    end
                end
            end
        end
    end



%% Smooth the streams aligned

for v=1:length(virus)    
    for d=1:length(dFF_names)
        for pow = 1:length(Opto_Powers)
            for k = 1:size(datatype,2)
                for nummice=1:size(fullstreams_aligned_pooled.(virus{v}).(Opto_Powers{pow}),1)
                                  
                    fullstreams_aligned_pooled_smooth.(virus{v}).(Opto_Powers{pow})(nummice,:) = smooth(fullstreams_aligned_pooled.(virus{v}).(Opto_Powers{pow})(nummice,:)',1000); %smooth 10 points at 100fps means smooth over 1 sec 
                
                end
                
            end
        end
    end
end

%% Downsample the streams aligned (smoothed)
for v=1:length(virus)    
    for d=1:length(dFF_names)
        for pow = 1:length(Opto_Powers)
            for k = 1:size(datatype,2)
                for nummice=1:size(fullstreams_aligned_pooled.(virus{v}).(Opto_Powers{pow}),1)
                
                    data2plot=fullstreams_aligned_pooled_smooth.(virus{v}).(Opto_Powers{pow})(nummice,:);
                    ttt = timetable(data2plot','SampleRate',sampling_rate_ds);
                    tttAvg = retime(ttt,'regular','mean','SampleRate',1); %1fps
                    fullstreams_aligned_pooled_ds.(virus{v}).(Opto_Powers{pow})(nummice,:) = tttAvg.Var1;
                end
                
                time_vect_aligned_ds.(Opto_Powers{pow}) = [seconds(tttAvg.Time)] + timevect_aligned.(virus{v}).(Opto_Powers{pow}).(AnimalIDs{nummice})(1); %also add the first timepoint from the previous timevector allowing to shift in time                
                time_vect_aligned_ds_minutes.(Opto_Powers{pow}) = time_vect_aligned_ds.(Opto_Powers{pow})/60-10;
                
            end
        end
    end
end



    
%% Plot 
color2plot = {'k','r','b','c','g','y','m'}; % food rest
% color2plot = {'r','b','k','c','g','y','m'}; % weeklz DETQ
figure; clf; hold on
if show_plot == 0
    set(gcf,'visible','off')
end

% AVERAGE GRAPH
idx_selct = [1,2,3,4,5,6,7]; % week1,2,3
% idx_selct = [1,2,3,4,5]; %fed vs food-restricted


for pow = 1:length(Opto_Powers)
    for v=1:length(virus)
        % average plot
        tmp_avg_pool = nanmean(fullstreams_aligned_pooled_ds.(virus{v}).(Opto_Powers{pow})(idx_selct,:),1); 
        tmp_error_pool = nanstd(fullstreams_aligned_pooled_ds.(virus{v}).(Opto_Powers{pow})(idx_selct,:),1,1)./...
        sqrt(sum(~isnan(fullstreams_aligned_pooled_ds.(virus{v}).(Opto_Powers{pow})(idx_selct,1))));
        error_area(time_vect_aligned_ds_minutes.(Opto_Powers{pow}),tmp_avg_pool,tmp_error_pool,color2plot{v},0.25); %t_trials is the time vector

    end
end

% Formatting
ylim([-10 60]); % food rest
% ylim([-20 80]); % weeklz DETQ

xlim([-10 50]);
xticks([-10 0 10 20 30 40 50]);
ax = gca;
ax.FontSize = 20; 
xlabel('Time (min)','FontSize', 20)
ylabel('\DeltaF/F (%)','FontSize', 20)
box off
set(gca,'TickLength',[0 0])

xline(-10,':k'); xline(-9,':k','Veh','FontSize',20); % grey out the noise in signal due to injection 
xline(0,':k'); xline(1,':k','DETQ','FontSize',20); % grey out the noise in signal due to injection 
% xline(20,':k'); xline(21,':k','J60 or Sal','FontSize',20); 
error_area_onlyrectangle([-10 -9],[40 120],[150 150],[165, 167, 168]./255,0.8); 
error_area_onlyrectangle([0 1],[40 120],[150 150],[165, 167, 168]./255,0.8); 
% error_area_onlyrectangle([20 21],[40 120],[150 150],[165, 167, 168]./255,0.8); 

%WEEKLY DETQ
yline(27.11,'-k'); yline(27.11*0.85,'--k'); yline(27.11*1.15,'--k','Plateau','FontSize',20,'FontAngle','italic'); % plateau
xline(5.238,'--k'); xline(25.1,'--k','FontSize',20); % temporal window
annotation('textbox', [0.348, 0.31, 0.2, 0.05], 'String', 'Temporal window', 'Color', 'black', 'FontSize', 18, 'EdgeColor', 'none','Interpreter', 'none','FontAngle','italic');

% week 1,2,3
annotation('textbox', [0.718, 0.90, 0.2, 0.05], 'String', 'Week 1', 'Color', 'red', 'FontSize', 18, 'EdgeColor', 'none','Interpreter', 'none');
annotation('textbox', [0.718, 0.85, 0.2, 0.05], 'String', 'Week 2', 'Color', 'blue', 'FontSize', 18, 'EdgeColor', 'none','Interpreter', 'none');
annotation('textbox', [0.718, 0.8, 0.2, 0.05], 'String', 'Week 3', 'Color', 'black', 'FontSize', 18, 'EdgeColor', 'none','Interpreter', 'none');

% % FED VS FOOD REST
% annotation('textbox', [0.60, 0.90, 0.2, 0.05], 'String', 'Ad lib fed', 'Color', 'black', 'FontSize', 18, 'EdgeColor', 'none','Interpreter', 'none','FitBoxToText','on');
% annotation('textbox', [0.60, 0.85, 0.2, 0.05], 'String', 'Food restriction', 'Color', 'red', 'FontSize', 18, 'EdgeColor', 'none','Interpreter', 'none','FitBoxToText','on');
% annotation('textbox', [0.60, 0.80, 0.2, 0.05], 'String', 'Ad lib fed', 'Color', 'blue', 'FontSize', 18, 'EdgeColor', 'none','Interpreter', 'none','FitBoxToText','on');

%saveplot or not
if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
    saveas(gcf,[PATH2DATA_0,'\pooled figures\streams_aligned_overlap ',' ',datatype{k},'.tif']);
    saveas(gcf,[PATH2DATA_0,'\pooled figures\streams_aligned_overlap ',' ',datatype{k},'.fig']);
    saveas(gcf,[PATH2DATA_0,'\pooled figures\streams_aligned_overlap ',' ',datatype{k},'.pdf']);
end          

%% SUBPLOT INDIDIVUAL GRAPH
color2plot = {'r','b','k','c','g','y','m'};

idx_selct = [1,2,3,4,5,6,7]; % week1,2,3
% idx_selct = [1,2,3,4,5]; %fed vs. food-restricted

figure; clf; hold on
if show_plot == 0
    set(gcf,'visible','off')
end

figure;

% Get screen size
screenSize = get(0, 'ScreenSize');

% Set the figure position to fill the entire screen
set(gcf, 'Position', [1, 1, screenSize(3)/2-20, screenSize(4)-80]);

for mouse_number = 1:length(idx_selct)
    subplot(4,2,mouse_number)
    for pow = 1:length(Opto_Powers)
        for v=length(virus):-1:1
            % average plot
            tmp_avg_pool = nanmean(fullstreams_aligned_pooled_ds.(virus{v}).(Opto_Powers{pow})(idx_selct(mouse_number),:),1); 
            tmp_error_pool = nanstd(fullstreams_aligned_pooled_ds.(virus{v}).(Opto_Powers{pow})(idx_selct(mouse_number),:),1,1)./...
            sqrt(sum(~isnan(fullstreams_aligned_pooled_ds.(virus{v}).(Opto_Powers{pow})(idx_selct(mouse_number),1))));
            error_area(time_vect_aligned_ds_minutes.(Opto_Powers{pow}),tmp_avg_pool,tmp_error_pool,color2plot{v},0.25); %t_trials is the time vector
        end
    end

% Formatting
ylim([-20 80]);  
xlim([-10 50]);
xticks([-10 0 10 20 30 40 50]);
ax = gca;
ax.FontSize = 16; 
xlabel('Time (min)','FontSize', 16)
ylabel('\DeltaF/F','FontSize', 16)
box off
set(gca,'TickLength',[0 0])

xline(-10,':k'); xline(-9,':k','V.','FontSize',16); % grey out the noise in signal due to injection 
xline(0,':k'); xl = xline(1,':k','D.','FontSize',16); % grey out the noise in signal due to injection 
% xl.LabelHorizontalAlignment = 'left';

% xline(20,':k'); xline(21,':k','J60 or Sal','FontSize',20); 
error_area_onlyrectangle([-10 -9],[40 120],[150 150],[165, 167, 168]./255,0.8); 
error_area_onlyrectangle([0 1],[40 120],[150 150],[165, 167, 168]./255,0.8); 
% error_area_onlyrectangle([20 21],[40 120],[150 150],[165, 167, 168]./255,0.8); 

mouse_number_string = string(mouse_number);
% Modify the title line to set the title centered above the subplot
t = title(['mouse ' char(mouse_number_string),':'], 'FontSize', 16, 'FontWeight', 'bold');

% Adjust the position of the title to be centered above the subplot
set(t, 'Position', [31, 68, 0]);

end

% Weekly DETQ
% Add text at the far bottom right of the figure with specified colors
annotation('textbox', [0.56, 0.25-0.015, 0.2, 0.05], 'String', 'Week 1', 'Color', 'red', 'FontSize', 16, 'EdgeColor', 'none');
annotation('textbox', [0.56, 0.225-0.015, 0.2, 0.05], 'String', 'Week 2', 'Color', 'blue', 'FontSize', 16, 'EdgeColor', 'none');
annotation('textbox', [0.56, 0.200-0.015, 0.2, 0.05], 'String', 'Week 3', 'Color', 'black', 'FontSize', 16, 'EdgeColor', 'none');

% Food restriction
% annotation('textbox', [0.56, 0.25-0.015, 0.2, 0.05], 'String', 'Fed', 'Color', 'red', 'FontSize', 16, 'EdgeColor', 'none');
% annotation('textbox', [0.56, 0.225-0.015, 0.2, 0.05], 'String', 'Food-restr', 'Color', 'blue', 'FontSize', 16, 'EdgeColor', 'none');
% annotation('textbox', [0.56, 0.200-0.015, 0.2, 0.05], 'String', 'Fed 2', 'Color', 'black', 'FontSize', 16, 'EdgeColor', 'none');

%saveplot or not
if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
    saveas(gcf,[PATH2DATA_0,'\pooled figures\indiv subplots ',' ',datatype{k},'.tif']);
    saveas(gcf,[PATH2DATA_0,'\pooled figures\indiv subplots ',' ',datatype{k},'.fig']);
    saveas(gcf,[PATH2DATA_0,'\pooled figures\indiv subplots ',' ',datatype{k},'.pdf']);

    % print([PATH2DATA_0,'\pooled figures\indiv subplots ',' ',datatype{k},'.pdf'], '-dpdf', '-fillpage');

end 


%% Quantifications

% idx_selct = [1,2,3,4,5,6,7]; %weekly DETQ
idx_selct = [1,2,3,4,5]; %fed vs. restricted

% getting the indexes
Section={'Drug1','Drug2'}
temp_t = time_vect_aligned_ds.VEH - time_vect_aligned_ds.VEH(1); 
dummie = 1:length(temp_t); 
dummie2 = dummie(temp_t >= 599); % 600 sec = 10 min --> first 10 min
idx_Drug1_2 = dummie2(1);  
idx_Drug1_1 = find(temp_t == 0);
idx_Drug1 = idx_Drug1_1:idx_Drug1_2; 

temp_t = time_vect_aligned_ds.DETQ - time_vect_aligned_ds.DETQ(1); 
dummie = 1:length(temp_t); 
dummie3 = dummie(temp_t >= 60*49-1);  % 49 min = end for weekly DETQ and for the plateau calc for food-restr; % 15 min for calculating the AUC food restriction
idx_Drug2_1 = dummie3(1); 
% dummie4 = find(temp_t >= 60); % 1min into DETQ
dummie4 = find(temp_t >= 0); % from time 0, in this case it's 1 min into DETQ (coded earlier in the extraction code)
idx_Drug2_2 = dummie4(1); 
idx_Drug2 = idx_Drug2_2:idx_Drug2_1-1; 

idx_drug_all = {'idx_Drug1','idx_Drug2'};
Day = {'veh','detq'};

% initialize Measurements
for du=1:length(Day)
    for v=1:length(virus)    
        for nummice2=1:length(idx_selct); %size(fullstreams_aligned_pooled.(virus{v}).(Opto_Powers{pow}),1)
            nummice = idx_selct(nummice2)

            Measurements_AVE.AUC.(Day{du}) = ones(length(virus),length(mice_list_virus.(virus{v})))*nan; %done like this, for easier copy-paste in GraphPad. day1 (rows=weeks (virus) and columns=animals)
  
        end                
    end
end

for du=1:length(Day) %here Day means drug (veh or detq)
    for v=1:length(virus) 
        for nummice2=1:length(idx_selct); %size(fullstreams_aligned_pooled.(virus{v}).(Opto_Powers{pow}),1)
            nummice = idx_selct(nummice2)
            
            if du == 1; idx_select = idx_Drug1;  pow = 1; elseif du == 2; idx_select = idx_Drug2; pow = 2; end
            % find the right dff using the indexes
            dff_select = fullstreams_aligned_pooled_ds.(virus{v}).(Opto_Powers{pow})(nummice,idx_select);

            % each row is a new week (here variable virus)
            % and each column is a new mouse
            Measurements_AVE.AUC.(Day{du})(v,nummice) = trapz(dff_select)/length(idx_select);

        end                
    end
end  
% end




%% Save pooled data for all viruses together
% STORE INDIVIDUAL and AVERAGE TRIALS
save([PATH2DATA_0,'\PooledAllMice.mat'],'virus','PooledAnimalID','dFF_names','datatype','Opto_Powers','fullstreams_aligned_pooled_smooth','Measurements_AVE',...
    'dt_ds','sampling_rate_ds','length_data','fullstreams','time_vect','pooledtype','fullstreams_aligned','timevect_aligned','fullstreams_aligned_pooled_ds','time_vect_aligned_ds','time_vect_aligned_ds_minutes');


    
    
  
        