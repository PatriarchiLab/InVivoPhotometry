%Averaging the data from the Unpredicted Reward experiment
%Author:Maria Wilhelm

%Clean the previous matlab sesion

close all
clearvars
clc

%give folder with matlab functions
folder_with_matlab_scripts = 'Y:\Chernysheva Maria\MATLAB\Matlab_scripts_AlloLight';
addpath (genpath(folder_with_matlab_scripts))

data_folder.Veh1 = 'Y:\Chernysheva Maria\Fiber_Photometry_data\PFC\DataForAveraging\Veh1_All_mice';
data_folder.DETQ = 'Y:\Chernysheva Maria\Fiber_Photometry_data\PFC\DataForAveraging\DETQ_All_mice';
data_folder.Veh2 = 'Y:\Chernysheva Maria\Fiber_Photometry_data\PFC\DataForAveraging\Veh2_All_mice';
Folder_4Figures = 'Y:\Chernysheva Maria\Fiber_Photometry_data\PFC\DataForAveraging';

RangeTrials = 12; %max ammoount of trials done by all mice in all sessions
for InjectionType = 1:length(fieldnames(data_folder)) %loop cohorts (Veh1, DETQ or Veh2 injections)
    
    folders = fieldnames(data_folder); %recover the name of the cohort from the Str
    cd  (data_folder.(folders{InjectionType}))% set the folder for extracting the data from the Cohort
    
    list_of_files = dir ('*.mat'); %list the data folder content
    
    for mouse_Number = 1:1:size(list_of_files,1) %loop mice
        name_file{mouse_Number} = list_of_files(mouse_Number).name; %Find the names of the files for each mouse
        load(name_file{mouse_Number}) %load the file
        data_injection_type{InjectionType}{mouse_Number} = StrStr;
        clearvars StrStr
        MouseID{mouse_Number} = fieldnames(data_injection_type{InjectionType}{mouse_Number});
        
        seconds_before = 5;
        seconds_after = 10;
        
        f_lowpass = 1;
        Trim = 60;
        dT_Doric_average = 0.0082; 
        FpS_perChannel = 1/dT_Doric_average;
        Time2Reward = [];
        Time2Reward = -1*(seconds_before):dT_Doric_average:seconds_after+dT_Doric_average;
        
        % Extract data
        index_Reward = data_injection_type{InjectionType}{mouse_Number}.(MouseID{mouse_Number}{1}).index_Reward;
        FunctionalSignal_demodulated_Clean = data_injection_type{InjectionType}{mouse_Number}.(MouseID{mouse_Number}{1}).FunctionalSignal_demodulated_Clean;
        IsosbesticSignal_demodulated_Clean = data_injection_type{InjectionType}{mouse_Number}.(MouseID{mouse_Number}{1}).IsosbesticSignal_demodulated_Clean;
        
        %calculate the number of data points incide the defined time
        index_Event = index_Reward;
        IndexBeforeEvent = round(seconds_before/dT_Doric_average);
        IndexAfterEvent = round(seconds_after/dT_Doric_average);
        dFF_G_filtered_Trial = [];
        Vector_IsosbesticSignal_Filtered_Trials = [];
        Vector_FunctionalSignal_Filtered_Trials = [];
        dFF_G_filtered_Trial_Trim =[];
        movmean_dFF_Trial_Trim = [];
        movmean_dFF_Trial_Trim_Norm = [];
        smooth_dFF_G_filtered_Trial_Trim = [];
        
        
        movWindow = 50;
        for ii=1:length(index_Event)
            
            Vector_FunctionalSignal_G_Trial = FunctionalSignal_demodulated_Clean(index_Event(ii)-(IndexBeforeEvent+Trim):index_Event(ii)+(IndexAfterEvent+Trim));
            Vector_IsosbesticSignal_G_Trial = IsosbesticSignal_demodulated_Clean(index_Event(ii)-(IndexBeforeEvent+Trim):index_Event(ii)+(IndexAfterEvent+Trim));
                        
            Vector_FunctionalSignal_G_Trial_Filtered = lowpass(Vector_FunctionalSignal_G_Trial,f_lowpass,FpS_perChannel);
            Vector_IsosbesticSignal_G_Trial_Filtered = lowpass(Vector_IsosbesticSignal_G_Trial,f_lowpass,FpS_perChannel);
            
            %fit the signal only during the baseline (before the event)
            X = Vector_IsosbesticSignal_G_Trial_Filtered(Trim+1:IndexBeforeEvent);
            Y = Vector_FunctionalSignal_G_Trial_Filtered(Trim+1:IndexBeforeEvent);
            
            PolynomialCoefficients_G = polyfit(X,Y,1);                     
            Y_fit_trial_G = PolynomialCoefficients_G(1) .* Vector_IsosbesticSignal_G_Trial_Filtered + PolynomialCoefficients_G(2);
            
            Vector_FunctionalSignal_G_Trial_Sub = Vector_FunctionalSignal_G_Trial_Filtered - Y_fit_trial_G;
            dFF_G_filtered_Trial(ii,:) = 100*((Vector_FunctionalSignal_G_Trial_Sub)./Y_fit_trial_G);
            %
            dFF_G_filtered_Trial_Trim(ii,:) = dFF_G_filtered_Trial(ii,Trim+1:length(dFF_G_filtered_Trial(ii,:))-Trim);
            smooth_dFF_G_filtered_Trial_Trim(ii,:) = movmean(dFF_G_filtered_Trial_Trim(ii,:),movWindow);
            Vector_IsosbesticSignal_Filtered_Trials(ii,:) = Vector_IsosbesticSignal_G_Trial_Filtered(Trim+1:length(dFF_G_filtered_Trial(ii,:))-Trim);
            Vector_FunctionalSignal_Filtered_Trials(ii,:) = Vector_FunctionalSignal_G_Trial_Filtered(Trim+1:length(dFF_G_filtered_Trial(ii,:))-Trim);
            
            Fbase_465 = mean(Vector_FunctionalSignal_G_Trial_Filtered(Trim+1:IndexBeforeEvent));
            Fbase_405 = mean(Vector_IsosbesticSignal_G_Trial_Filtered(Trim+1:IndexBeforeEvent));
            
            dFF_G_filtered_Trial_465(ii,:) = 100*(Vector_FunctionalSignal_G_Trial_Filtered'-Fbase_465)/Fbase_465;
            dFF_G_filtered_Trial_405(ii,:) = 100*(Vector_IsosbesticSignal_G_Trial_Filtered'-Fbase_405)/Fbase_405;
            
            dFF_G_filtered_Trial_Trim_465(ii,:) = dFF_G_filtered_Trial_465(ii,Trim+1:length(dFF_G_filtered_Trial(ii,:))-Trim);
            dFF_G_filtered_Trial_Trim_405(ii,:) = dFF_G_filtered_Trial_405(ii,Trim+1:length(dFF_G_filtered_Trial(ii,:))-Trim);
            
            z_score_dFF(ii,:) = zscore_norm(dFF_G_filtered_Trial_Trim(ii,:), 1, IndexBeforeEvent);
            
            z_score_dFF_smoothed(ii,:) = movmean (z_score_dFF(ii,:),movWindow);
           
        end
        dFF_G_filtered_Trial_Trim_allMice{InjectionType}{mouse_Number} = dFF_G_filtered_Trial_Trim;

        z_score_filtered_Trial_smoothed_dFF{InjectionType}{mouse_Number} = z_score_dFF_smoothed;
        
        
        dFF_G_filtered_Trial_Trim_465_allMice{InjectionType}{mouse_Number} = dFF_G_filtered_Trial_Trim_465;
        dFF_G_filtered_Trial_Trim_405_allMice{InjectionType}{mouse_Number} = dFF_G_filtered_Trial_Trim_405;
        
        Isosbestic_filtered_Trial_Trim_allMice{InjectionType}{mouse_Number} = Vector_IsosbesticSignal_Filtered_Trials;
        Functional_filtered_Trial_Trim_allMice{InjectionType}{mouse_Number} = Vector_FunctionalSignal_Filtered_Trials;
        
        
        smooth_dFF_G_filtered_Trial_Trim_allMice{InjectionType}{mouse_Number} = smooth_dFF_G_filtered_Trial_Trim;        
        dFF_G_filtered_Trial_Trim_MeanPerMouse{InjectionType}{mouse_Number} = mean(dFF_G_filtered_Trial_Trim_allMice{InjectionType}{mouse_Number}(1:RangeTrials,:),1);        
        dFF_G_filtered_Trial_Trim_MeanPerMouse_465{InjectionType}{mouse_Number} = mean(dFF_G_filtered_Trial_Trim_465_allMice{InjectionType}{mouse_Number}(1:RangeTrials,:),1);
        dFF_G_filtered_Trial_Trim_MeanPerMouse_405{InjectionType}{mouse_Number} = mean(dFF_G_filtered_Trial_Trim_405_allMice{InjectionType}{mouse_Number}(1:RangeTrials,:),1);      
        
        dFF_G_filtered_Trial_Trim_MeanPerMouse_movmean{InjectionType}{mouse_Number} = movmean(dFF_G_filtered_Trial_Trim_MeanPerMouse{InjectionType}{mouse_Number},movWindow);        
        smooth_dFF_G_filtered_Trial_Trim_MeanPerMouse{InjectionType}{mouse_Number} = mean(smooth_dFF_G_filtered_Trial_Trim_allMice{InjectionType}{mouse_Number}(1:RangeTrials,:),1);              
        Peak_Reward{InjectionType}{mouse_Number} = max(dFF_G_filtered_Trial_Trim_MeanPerMouse_movmean{InjectionType}{mouse_Number}(IndexBeforeEvent:IndexBeforeEvent+round(2/dT_Doric_average)));
        Peak_Reward_relative_to_noise{InjectionType}{mouse_Number} = mean(dFF_G_filtered_Trial_Trim_MeanPerMouse_movmean{InjectionType}{mouse_Number}(IndexBeforeEvent:IndexBeforeEvent+round(2/dT_Doric_average)))/mean(dFF_G_filtered_Trial_Trim_MeanPerMouse_movmean{InjectionType}{mouse_Number}(IndexBeforeEvent-round(4/dT_Doric_average):IndexBeforeEvent-round(2/dT_Doric_average)));
    end
    
end

%% find mean duration of the reward consumption
d_Index_Rewarded_Licking=[];
d_Index_Rewarded_Licking_Start=[];
d_Index_Rewarded_Licking_Finish=[];
for InjectionType=1:length(fieldnames(data_folder))
    for mouse_Number = 1:1:size(list_of_files,1)
        MouseID{mouse_Number} = fieldnames(data_injection_type{InjectionType}{mouse_Number});
        
        index_Reward_all{InjectionType}{mouse_Number} = data_injection_type{InjectionType}{mouse_Number}.(MouseID{mouse_Number}{1}).index_Reward;
        index_startlick_all{InjectionType}{mouse_Number} = data_injection_type{InjectionType}{mouse_Number}.(MouseID{mouse_Number}{1}).index_startlick;
        index_finishlick_all{InjectionType}{mouse_Number} = data_injection_type{InjectionType}{mouse_Number}.(MouseID{mouse_Number}{1}).index_finishlick;
        licking_vector{InjectionType}{mouse_Number} = data_injection_type{InjectionType}{mouse_Number}.(MouseID{mouse_Number}{1}).LickingResponse_Clean;
        
        total_licking_percent{InjectionType}{mouse_Number} = 100*(length(nonzeros(licking_vector{InjectionType}{mouse_Number}))/length(licking_vector{InjectionType}{mouse_Number}));
        [d_Index_Rewarded_Licking{InjectionType}{mouse_Number},d_Index_Rewarded_Licking_Start{InjectionType}{mouse_Number},d_Index_Rewarded_Licking_Finish{InjectionType}{mouse_Number}]=LickingDuration(index_Reward_all{InjectionType}{mouse_Number}, index_finishlick_all{InjectionType}{mouse_Number},index_startlick_all{InjectionType}{mouse_Number});
    end
end

d_Index_Rewarded_Licking_allcells=[d_Index_Rewarded_Licking{:}];
d_Index_Rewarded_Licking_Start_allcells=[d_Index_Rewarded_Licking_Start{:}];
d_Index_Rewarded_Licking_Finish_allcells=[d_Index_Rewarded_Licking_Finish{:}];

mean_dT_Rewarded_Licking=mean([d_Index_Rewarded_Licking_allcells{:}])*dT_Doric_average;
mean_dT_Rewarded_Licking_Start=mean([d_Index_Rewarded_Licking_Start_allcells{:}])*dT_Doric_average;
mean_dT_Rewarded_Licking_Finish=mean([d_Index_Rewarded_Licking_Finish_allcells{:}])*dT_Doric_average;
%% quantify licking response
Color_Grey = [0.7 0.7 0.7];
Color_NiceRed = [255, 128, 128]./255;
Color_Veh1 = [0.2 0.2 0.2];
Color_DETQ = Color_NiceRed;
Color_Veh2 = Color_Grey;

name_plot = 'Total Licking Duration';
fontSize_4figures=12;
input1 = round(vertcat(total_licking_percent{1}{:}),3);
input2 = round(vertcat(total_licking_percent{2}{:}),3);
input3 = round(vertcat(total_licking_percent{3}{:}),3);

paired=1;
figure
bar_plots_open_field({input1,input2,input3},{LegendVeh1;LegendDETQ;LegendVeh},{Color_Veh1;Color_DETQ;Color_Veh2}, name_plot,'Licking Duration (%)',fontSize_4figures, 1)
save_plots(name_plot)

excel_input = horzcat(input1,input2,input3);
%%
NMice=size(list_of_files,1);

LineWidth=1.2;

LegendVeh1=['\color[rgb]{',num2str(Color_Veh1),'}Veh'];
LegendDETQ=['\color[rgb]{',num2str(Color_DETQ),'}DETQ'];
LegendVeh=['\color[rgb]{',num2str(Color_Veh2),'}Veh'];

Color_Lick=[0.6431    0.4549    0.2863 0.5];
Color_Lick=Color_Lick(1:3);
NamePlot=['Unexpected Reward _','z-score','12 Trials per Mouse','5-10 seconds','movmean50','Design2'];

Input2Plot=dFF_G_filtered_Trial_Trim_MeanPerMouse;
%
input_1=vertcat(Input2Plot{1}{:});
input_2=vertcat(Input2Plot{2}{:});
input_3=vertcat(Input2Plot{3}{:});

Mean_deltaFF_Day1=mean(input_1,1);
sd_Mean_deltaFF_Day1=std(input_1,1)/sqrt(size(input_1,1));

Mean_deltaFF_Day2=mean(input_2,1);
sd_Mean_deltaFF_Day2=std(input_2,1)/sqrt(size(input_2,1));

Mean_deltaFF_Day3=mean(input_3,1);
sd_Mean_deltaFF_Day3=std(input_3,1)/sqrt(size(input_3,1));

figure

hold on

H1 = shadedErrorBar(Time2Reward,movmean(Mean_deltaFF_Day1,movWindow),movmean(sd_Mean_deltaFF_Day1,movWindow),{'color', Color_Veh1,'linewidth',LineWidth},0);
hold on
H2 = shadedErrorBar(Time2Reward,movmean(Mean_deltaFF_Day2,movWindow),movmean(sd_Mean_deltaFF_Day2,movWindow),{'color', Color_DETQ,'linewidth',LineWidth},0);
hold on
H3 = shadedErrorBar(Time2Reward,movmean(Mean_deltaFF_Day3,movWindow),movmean(sd_Mean_deltaFF_Day3,movWindow),{'color', Color_Veh2,'linewidth',LineWidth},0);
main_lines{1} = H1.mainLine;
main_lines{2} = H1.mainLine;
main_lines{3} = H1.mainLine;

hold on
fillearea=fill ([-1*(mean_dT_Rewarded_Licking_Start) 1*(mean_dT_Rewarded_Licking_Finish) 1*(mean_dT_Rewarded_Licking_Finish) -1*(mean_dT_Rewarded_Licking_Start)],[-0.8 -0.8 -1 -1],Color_Lick,'edgecolor', 'none');
fillearea.FaceAlpha = 0.5;

title (['Unexpected Reward','',{''}],'Interpreter', 'none')
ylabel('\DeltaF/F (%)','Fontsize',14,'Color','k');
% ylabel('z-score','Fontsize',14,'Color','k');
xlabel('Time to Reward, s','Fontsize',14,'Color','k');
hold on
xline(0,':k')

legend({LegendVeh1, ['\color[rgb]{',num2str(Color_DETQ),'}DETQ'], ['\color[rgb]{',num2str(Color_Veh2),'}Veh'], 'Licking'},'Location','eastoutside','FontSize',16)
set(gca,'FontName','Arial');
set(gca,'FontSize',16);
set(gca,'TickDir','out');
set(gca,'TickLength',[0.03, 0.1])
pbaspect([1.5 1 1])
box off
legend boxoff

cd (Folder_4Figures)
print(gcf, NamePlot,'-dpdf')
savefig([NamePlot,'.fig'])
%% find proportion of trials where the peak upon the reward is detected
input_for_single_trial_analysis = z_score_filtered_Trial_smoothed_dFF;

for InjectionType = 1:size (input_for_single_trial_analysis,2)
    
    for current_mouse = 1:size (input_for_single_trial_analysis{1, InjectionType},2)
        for current_trial = 1:size(input_for_single_trial_analysis{InjectionType}{current_mouse},1)
            [peak_prominence{InjectionType}{current_mouse}{current_trial},...
                peak_coordinate{InjectionType}{current_mouse}{current_trial},...
                w,p] = findpeaks (input_for_single_trial_analysis{InjectionType}{current_mouse}(current_trial,(IndexBeforeEvent:IndexBeforeEvent+round(2/dT_Doric_average))),'MinPeakProminence',1);
            
            [peak_prominence_baseline{InjectionType}{current_mouse}{current_trial},...
                peak_baseline{InjectionType}{current_mouse}{current_trial},...
                w_baseline,p_baseline] = findpeaks (input_for_single_trial_analysis{InjectionType}{current_mouse}(current_trial,(IndexBeforeEvent-round(2/dT_Doric_average):IndexBeforeEvent)),'MinPeakProminence',1);
        
        detected_peaks_per_trial{InjectionType}{current_mouse} = find(~cellfun(@isempty,peak_coordinate{InjectionType}{current_mouse}));
        detected_peaks_per_trial_baseline{InjectionType}{current_mouse} = find(~cellfun(@isempty,peak_coordinate{InjectionType}{current_mouse}));
        normalized_peak_prominence{InjectionType}{current_mouse}(current_trial) = mean(peak_prominence{InjectionType}{current_mouse}{current_trial})/mean(peak_prominence_baseline{InjectionType}{current_mouse}{current_trial});
        end
        normalized_peak_prominence_per_mouse{InjectionType}(current_mouse) = mean (normalized_peak_prominence{InjectionType}{current_mouse}(~isnan(normalized_peak_prominence{InjectionType}{current_mouse})));
        number_of_peaks_per_session {InjectionType}{current_mouse} = length(detected_peaks_per_trial{InjectionType}{current_mouse});
        percent_of_detected_peaks_per_session {InjectionType}(current_mouse)= 100*(number_of_peaks_per_session {InjectionType}{current_mouse}/size(input_for_single_trial_analysis{InjectionType}{current_mouse},1));
    end
end
NamePlot='Percent of Reward Response';
FontSize=12;


input1 = (percent_of_detected_peaks_per_session {1})';
input2 = (percent_of_detected_peaks_per_session {2})';
input3 = (percent_of_detected_peaks_per_session {3})';
paired=1;
figure
SEM_Input = bar_plots_open_field({input1,input2,input3},{LegendVeh1;LegendDETQ;LegendVeh},{Color_Veh1;Color_DETQ;Color_Veh2}, NamePlot,'Peaks to Reward (%)',FontSize, paired);
print(gcf, NamePlot,'-dpdf')
savefig([NamePlot,'.fig'])

excel_source_data_fig6_panel_j = vertcat ([input1, input2, input3],...
                                              [mean(input1),mean(input2),mean(input3)],...
                                               SEM_Input);

%%
%source data
for current_group = 1:size (dFF_G_filtered_Trial_Trim_MeanPerMouse,2)
    
    for current_mouse = 1:size (dFF_G_filtered_Trial_Trim_MeanPerMouse{1, current_group},2)
        
        dFF_MeanPerMouse_for_source_data {current_group}{current_mouse} = dFF_G_filtered_Trial_Trim_MeanPerMouse{1, current_group}{1, current_mouse}';
        
    end
end

Time2Reward_for_source_data = Time2Reward';

Mean_deltaFF_Day1_for_source_data = movmean(Mean_deltaFF_Day1',movWindow);
SEM_deltaFF_Day1_for_source_data = movmean(sd_Mean_deltaFF_Day1',movWindow);

Mean_deltaFF_Day2_for_source_data = movmean(Mean_deltaFF_Day2',movWindow);
SEM_deltaFF_Day2_for_source_data = movmean(sd_Mean_deltaFF_Day2',movWindow);

Mean_deltaFF_Day3_for_source_data = movmean(Mean_deltaFF_Day3',movWindow);
SEM_deltaFF_Day3_for_source_data = movmean(sd_Mean_deltaFF_Day3',movWindow);

%% calculate peak dFF
NamePlot='Peak Fluorescence 12 trials 50 movemean';
FontSize=12;
input1=round(vertcat(Peak_Reward{1}{:}),3);
input2=round(vertcat(Peak_Reward{2}{:}),3);
input3=round(vertcat(Peak_Reward{3}{:}),3);
[p_value,h]=ranksum(input1,input2);
paired=1;
figure
bar_plots_open_field({input1,input2,input3},{LegendVeh1;LegendDETQ;LegendVeh},{Color_Veh1;Color_DETQ;Color_Veh2}, NamePlot,'\DeltaF/F (%)',FontSize, paired)

savefig([NamePlot,'.fig'])
excel_input = horzcat(input1,input2,input3);

%% plot the heat maps

NamePlot='Heat Map 06 Red';
clims=[0 0.6];

F=figure;

% C=dFF_G_filtered_Trial_Trim_allMice;
C=smooth_dFF_G_filtered_Trial_Trim_allMice;
x=Time2Reward;
y=(1:size(C,1));

subplot(1,3,1)

imagesc(x,y,C{1}{1}(1:RangeTrials,:),clims);
xlabel('Time to Reward (s)', 'FontSize', 16);
hold on
xline(0,'w','Linewidth',2)
title(['Veh',{''}], 'FontSize', 16,'Color',Color_Veh1);
set(gca,'FontName','Arial');
set(gca,'FontSize',16);
ylabel('Trials', 'FontSize', 20);
set(gca,'TickDir','out');
set(gca,'TickLength',[0.02, 0.1])
box off

hold on
s2=subplot(1,3,2);
imagesc(x,y,C{2}{1}(1:RangeTrials,:),clims);
box off
hold on
xline(0,'w','Linewidth',2)
title(['DETQ',{''}],'Color',Color_DETQ);
xlabel('Time to Reward (s)', 'FontSize', 16);
set(gca,'FontName','Arial');
% set(gca,'FontSize',16);
set(gca,'TickDir','out');
set(gca,'TickLength',[0.02, 0.1])
box off
set(gca,'Yticklabel',[])

hold on
s3=subplot(1,3,3);
imagesc(x,y,C{3}{1}(1:RangeTrials,:),clims);
hold on
xline(0,'w','Linewidth',2)
title(['Veh',{''}],'Color',Color_Veh2);
xlabel('Time to Reward (s)', 'FontSize', 16);
box off
set(gca,'Yticklabel',[])

set(gca,'FontName','Arial');
set(gca,'FontSize',16);
set(gca,'TickDir','out');
set(gca,'TickLength',[0.02, 0.1])

h=colorbar;
h_ylablel=ylabel(h,'\DeltaF/F (%)','Rotation',270.0,'FontSize',20)
h.Label.Position(1) = 4;
set(gcf, 'Position', get(0, 'Screensize'))
set(gcf,'PaperOrientation','landscape');

set(s2,'Position',[0.3808    0.1100    0.2134    0.8150],'FontSize',16);

set(s3,'Position',[0.6316    0.1100    0.2134    0.8150],'FontSize',16);
%
savefig([NamePlot,'.fig'])


close all
excel_source_data_fig6k = horzcat(Time2Reward',...
                                    nan(length(Time2Reward),1),...
                                    C{1}{1}(1:RangeTrials,:)',...
                                         nan(length(Time2Reward),1),...
                                  C{2}{1}(1:RangeTrials,:)',...
                                         nan(length(Time2Reward),1),...
                                  C{3}{1}(1:RangeTrials,:)');
