%averaging Open field data. Similar to dLight, just arrange to pic up data
%for AlloLightCtrl

%this function allows to read photometry data collected during Open Field
%Test, detrend and calculate dFF, as well as detect peack and quantify
%AUC

%Author: Maria Wilhelm

%Clean the previous matlab sesion:
close all
clearvars
clc

%give folder with matlab functions
folder_with_matlab_scripts = 'Y:\Chernysheva Maria\MATLAB\Matlab_scripts_AlloLight';
addpath (genpath(folder_with_matlab_scripts))

%give the folder with the data subfolders
main_folder_with_data_subfolders = 'Y:\Chernysheva Maria\Fiber_Photometry_data\data_for_averaging\Open_Field\AlloLite_Ctr';
cd (main_folder_with_data_subfolders)

set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 5 20 10]);% define the default size for the matlab figure
%Average the mean trace profile (B)
group_name{1} = 'veh_day1'; 
group_name{2} = 'detq_day2'; 
group_name{3} = 'veh_day3'; 

%indicate minutes from the start of the recording (recorded in the experimental log
%book)
time_ip_veh_day1=[33,  25,    27,   26]; %m008, m012, m018, m020, m021, m022
time_ip_detq_day2=[33,  25,    27,   26]; %m008, m012, m018, m020, m021, m022
time_ip_veh_day3=[33,  25,    27,   26]; %m008, m012, m018, m020, m021, m022


index_ip.veh_day1 = round(time_ip_veh_day1*60/0.0082);
index_ip.detq_day2 = round(time_ip_detq_day2*60/0.0082);
index_ip.veh_day3 = round(time_ip_veh_day3*60/0.0082);

index_100minutes = round(100*60/0.0082);
index_20minutes = round(20*60/0.0082);
index_22minutes = round(22*60/0.0082);
index_90minutes = round(90*60/0.0082);

%% 
number_of_different_groups = 3; 

for group_number=1:number_of_different_groups
folders_with_data{group_number} = uigetdir(main_folder_with_data_subfolders,['Select Folder with the Drug ',num2str(group_number)]);
end

%set the folder where you would like to store Figures and Variables for
%this experiment: 

folder_4_figures = uigetdir(main_folder_with_data_subfolders,'Select Folder to save Plots and Vars');

cd (folder_4_figures)
%% 
for group_number = 1:number_of_different_groups %loop drugs (saline vs yohimbine etc)
    
cd  (folders_with_data{group_number})% set the folder for extracting the data from the Cohort

list_of_files = dir ('*.mat'); %list the data folder content

% Make a structure with all variables
number_of_mice{group_number} = size(list_of_files,1);

for mouse_number = 1:1:number_of_mice{group_number}
   %Find the names of the files for each mouse 
   file_name{group_number}{mouse_number} = list_of_files(mouse_number).name(1:end-4); 

   data_cohort{group_number}{mouse_number} = load(file_name{group_number}{mouse_number}) ;
   
end
end

%% 
n_polyfit = 3;
window_median = 1; 

for group_number = 1:number_of_different_groups

for mouse_number = 1:1:number_of_mice{group_number}

current_mouse_data = data_cohort{1,group_number}{1,mouse_number}.(char(fieldnames(data_cohort{1,group_number}{1,mouse_number})));
current_index = index_ip.(group_name{group_number});
current_mouse_ID = char(fieldnames(data_cohort{1,group_number}{1,mouse_number}));
current_mouse_ID = current_mouse_ID(1:4);
mice_IDs {group_number}{mouse_number} = current_mouse_ID;
current_functional_signal = movmedian (current_mouse_data.FunctionalSignal_demodulated_Clean  , window_median);
current_isosbestic_signal = movmedian (current_mouse_data.IsosbesticSignal_demodulated_Clean, window_median);
current_time_vector = current_mouse_data.ConsoleTime_Clean  ;

index_end_fit = length (current_functional_signal);
index_start_fit = current_index(mouse_number)-index_22minutes;
index_range_baseline = horzcat([index_start_fit:1:current_index(mouse_number)],[current_index(mouse_number)+index_90minutes-index_20minutes:1:index_end_fit]);

baseline_signal_2fit_functional = current_functional_signal(index_range_baseline);

fit_2plot_functional = polyfit(current_time_vector (index_range_baseline),baseline_signal_2fit_functional,n_polyfit);

current_y_fit_functional = polyval(fit_2plot_functional,current_time_vector);
functional_signal_detrended = current_functional_signal - current_y_fit_functional;

baseline_signal_2fit_isosbestic = current_isosbestic_signal(index_range_baseline);
fit_2plot_isosbestic = polyfit(current_time_vector (index_range_baseline),baseline_signal_2fit_isosbestic,n_polyfit);
current_y_fit_isosbestic = polyval(fit_2plot_isosbestic,current_time_vector);

isosbestic_signal_detrended = current_isosbestic_signal - current_y_fit_isosbestic;
dFF_functional = current_functional_signal./current_y_fit_functional;
dFF_isosbestic = current_isosbestic_signal./current_y_fit_isosbestic;

dFF_normalized  = 100*((dFF_functional)-(dFF_isosbestic));

index_range_2plot = [current_index(mouse_number) - index_20minutes:1:current_index(mouse_number)+index_90minutes];

dFF_normalized_trimmed {1,group_number}{1,mouse_number} = dFF_normalized(index_range_2plot).';
smooth_dFF_normalized_trimmed {1,group_number}{1,mouse_number} = movmean(dFF_normalized_trimmed {1,group_number}{1,mouse_number},50);

%use moving median filter to calculate the baseline signal level
movmedian_window = 10000;
dFF_normalized_trimmed_baseline {1,group_number}{1,mouse_number} = movmedian(dFF_normalized_trimmed {1,group_number}{1,mouse_number},movmedian_window);

functional_signal.(group_name{group_number}).(mice_IDs {group_number}{mouse_number}) = current_functional_signal;
y_fit_functional.(group_name{group_number}).(mice_IDs {group_number}{mouse_number}) = current_y_fit_functional;
isosbestic_signal.(group_name{group_number}).(mice_IDs {group_number}{mouse_number}) = current_isosbestic_signal;
y_fit_isosbestic.(group_name{group_number}).(mice_IDs {group_number}{mouse_number}) = current_y_fit_isosbestic;
    end
end
%% plot average
time_trial=-20*60:0.0082:0.0082*index_90minutes;
Color_LightBlue=[ 0.5843    0.8157    0.9882];
figureColor_LightBlue=[ 0.5843    0.8157    0.9882];
MethodSmooth = 'movmean';

line_width = 1;
transparent = 1;
Color_Grey=[0.7 0.7 0.7];clc

Color_NiceRed=[255, 128, 128]./255;

Color_Veh1=[0.2 0.2 0.2];
Color_DETQ=Color_NiceRed;
Color_Veh2=Color_Grey;

Color_Plot{1} = Color_Veh1;
Color_Plot{2} = Color_DETQ;
Color_Plot{3} = Color_Veh2;

dFF_Input_mean_mouse = cell(1,number_of_different_groups);
dFF_G_filtered_AllMice_Averages = cell(1,number_of_different_groups);
Mean_dFF_AverageInput = cell(1,number_of_different_groups);
sd_Mean_dFF_AverageInput = cell(1,number_of_different_groups);
for group_number=1:number_of_different_groups
    
    for mouse_number = 1:1:number_of_mice{group_number}
    %average
    
    dFF_Input_mean_mouse{group_number}{mouse_number} = dFF_normalized_trimmed_baseline{1,group_number}{1,mouse_number};
    
    end
    
    dFF_G_filtered_AllMice_Averages{group_number}=vertcat(dFF_Input_mean_mouse{group_number}{:});   
    
    Mean_dFF_AverageInput{group_number} = mean(dFF_G_filtered_AllMice_Averages{group_number},1);
    sd_Mean_dFF_AverageInput{group_number}=std(dFF_G_filtered_AllMice_Averages{group_number},1)/sqrt(size(dFF_G_filtered_AllMice_Averages{group_number},1));    
   
end
%% 
step_2plot = 400; %to reduce the output pdf file size let's plot each50th point
Window_smooth = step_2plot;
figure

for group_number = 1:number_of_different_groups
    
mean_input_2plot = smoothdata(Mean_dFF_AverageInput{group_number},MethodSmooth,Window_smooth);
sd_input_2plot = smoothdata(sd_Mean_dFF_AverageInput{group_number},MethodSmooth,Window_smooth);
H = shadedErrorBar(time_trial(1:step_2plot:end)/60,mean_input_2plot(1:step_2plot:end),sd_input_2plot(1:step_2plot:end),{'color', Color_Plot{group_number},'linewidth',line_width},transparent);
hold on
end
xline(0,':')
hold off

fontSize_4figures=12;
name_plot='Average Baseline 90 min  after ip Min';
title ([name_plot,{''}])
ylabel('\DeltaF/F_0 (%)','Fontsize',24,'Color','k');
xlabel('Time to injection, min','Fontsize',fontSize_4figures,'Color','k');

xticks([-20:10:90])

xlim([-20 90])
ylim ([-5 25])
box off
set(gca,'TickDir','out');
set(gca,'FontSize',fontSize_4figures)
set(gca, 'Color', 'None');

cd (folder_4_figures)
save_plots(name_plot)

%% heat maps
NamePlot = 'Heat Map 40 Red';
color_limits = [0 8];

F = figure;

for group_number = 1:number_of_different_groups
% C=dFF_G_filtered_Trial_Trim_allMice;
C = dFF_G_filtered_AllMice_Averages{group_number};
x = time_trial/60;
y = (1:size(C,1));

subplot(1,number_of_different_groups,group_number)

imagesc(x,y,C(1:number_of_mice{group_number},:),color_limits);
xlabel('Time to ip (s)', 'FontSize', 14);
hold on
colorbar ('Ticks',[color_limits(1):1:color_limits(2)])

xline(0,'w','Linewidth',2)
title([group_name{group_number},{''}], 'FontSize', 14,'Color',Color_Plot{group_number}, 'Interpreter', 'none');
set(gca,'FontName','Arial');
set(gca,'FontSize',16);
ylabel('Mouse #', 'FontSize', 14);
set(gca,'TickDir','out');
set(gca,'TickLength',[0.02, 0.1])
box off

end

set (gcf, 'Position', [5.0000    5.0000   20.0000    7.4619]);
name_plot = 'Heat Maps abs axis';
save_plots(name_plot)

%% calculate AUC
cd (folder_4_figures)
AUC_start_time_after_ip_minutes = 5;
AUC_finish_time_after_ip_minutes = 40;

name_plot = ['AUC of Baseline ', num2str(AUC_start_time_after_ip_minutes),'-',num2str(AUC_finish_time_after_ip_minutes),' min'];

index_ip = index_20minutes;
index_start_AUC = index_ip + round(AUC_start_time_after_ip_minutes*60/0.0082);
index_finish_AUC = index_ip + round(AUC_finish_time_after_ip_minutes*60/0.0082);

for group_number = 1:number_of_different_groups
    
    for mouse_number = 1:1:number_of_mice{group_number}
    
        AUC_dFF{group_number}{mouse_number} = trapz(time_trial(index_start_AUC:index_finish_AUC)/60,dFF_normalized_trimmed_baseline{1,group_number}{1,mouse_number}(index_start_AUC:index_finish_AUC));
    
    end
    
    AUC_AllMice_Averages{group_number}=vertcat(AUC_dFF{group_number}{:});   
end
 


LegendVeh1=['\color[rgb]{',num2str(Color_Veh1),'}Veh'];
LegendDETQ=['\color[rgb]{',num2str(Color_DETQ),'}DETQ'];
LegendVeh=['\color[rgb]{',num2str(Color_Veh2),'}Veh'];


input1 = AUC_AllMice_Averages{1};
input2 = AUC_AllMice_Averages{2};
input3 = AUC_AllMice_Averages{3};

figure
bar_plots_open_field({input1,input2,input3},{LegendVeh1;LegendDETQ;LegendVeh},{Color_Veh1;Color_DETQ;Color_Veh2}, name_plot,'AUC, a.u.',fontSize_4figures, 1)
ylim ([-100 800]);

save_plots(name_plot)

%% find peaks 
peak_prominence_thresholds = 0.5;
signal_for_peak_detection = smooth_dFF_normalized_trimmed;

AUC_start_time_after_ip_minutes = 5;
AUC_finish_time_after_ip_minutes = 40;
index_ip_average = index_20minutes;
index_start_AUC = index_ip_average + round(AUC_start_time_after_ip_minutes*60/0.0082);
index_finish_AUC = index_ip_average + round(AUC_finish_time_after_ip_minutes*60/0.0082);
range_to_quantify = (index_start_AUC:1:index_finish_AUC);

for group_number = 1:number_of_different_groups
for mouse_number = 1:1:number_of_mice{group_number}
[peak_local_maxima.after_ip_5_40_minutes.(['peak_prominence_',num2str(round(peak_prominence_thresholds))]){1,group_number}{1,mouse_number},...
 peak_location.after_ip_5_40_minutes.(['peak_prominence_',num2str(round(peak_prominence_thresholds))]){1,group_number}{1,mouse_number},...
 w,...
 peak_prominence.after_ip_5_40_minutes.(['peak_prominence_',num2str(round(peak_prominence_thresholds))]){1,group_number}{1,mouse_number}] = findpeaks(signal_for_peak_detection{1,group_number}{1,mouse_number}(range_to_quantify),time_trial(range_to_quantify)/60,'MinPeakProminence',peak_prominence_thresholds);


[peak_local_maxima.baseline.(['peak_prominence_',num2str(round(peak_prominence_thresholds))]){1,group_number}{1,mouse_number},...
 peak_location.baseline.(['peak_prominence_',num2str(round(peak_prominence_thresholds))]){1,group_number}{1,mouse_number},...
 w,...
 peak_prominence.baseline.(['peak_prominence_',num2str(round(peak_prominence_thresholds))]){1,group_number}{1,mouse_number}] = findpeaks(signal_for_peak_detection{1,group_number}{1,mouse_number}(1:index_20minutes),time_trial(1:index_20minutes)/60,'MinPeakProminence',peak_prominence_thresholds);

number_of_peaks.after_ip_5_40_minutes(mouse_number,group_number) = length(peak_prominence.after_ip_5_40_minutes.(['peak_prominence_',num2str(round(peak_prominence_thresholds))]){1, group_number}{1, mouse_number}  );
number_of_peaks.baseline(mouse_number,group_number) = length(peak_prominence.baseline.(['peak_prominence_',num2str(round(peak_prominence_thresholds))]){1, group_number}{1, mouse_number}  );

amplitude_of_peaks.after_ip_5_40_minutes(mouse_number,group_number) = mean(peak_prominence.after_ip_5_40_minutes.(['peak_prominence_',num2str(round(peak_prominence_thresholds))]){1, group_number}{1, mouse_number}  );
amplitude_of_peaks.baseline(mouse_number,group_number) = mean(peak_prominence.baseline.(['peak_prominence_',num2str(round(peak_prominence_thresholds))]){1, group_number}{1, mouse_number}  );
end


end

% make different versions of the plots
%barplots with #N peaks and mean Amplitude

cd (folder_4_figures)
name_plot = ['Number of ', ' dFF Peaks'];

LegendVeh1=['\color[rgb]{',num2str(Color_Veh1),'}Veh'];
LegendDETQ=['\color[rgb]{',num2str(Color_DETQ),'}DETQ'];
LegendVeh=['\color[rgb]{',num2str(Color_Veh2),'}Veh'];
%% 

% % 
input1 = number_of_peaks.after_ip_5_40_minutes (:,1);
input2 = number_of_peaks.after_ip_5_40_minutes (:,2);
input3 = number_of_peaks.after_ip_5_40_minutes (:,3);


figure
bar_plots_open_field({input1,input2,input3},{LegendVeh1;LegendDETQ;LegendVeh},{Color_Veh1;Color_DETQ;Color_Veh2}, name_plot,'Number of Peaks',fontSize_4figures, 1)
save_plots(name_plot)

%% 
cd (folder_4_figures)
name_plot = 'Average Amplitude';

input1 = amplitude_of_peaks.after_ip_5_40_minutes (:,1);
input2 = amplitude_of_peaks.after_ip_5_40_minutes (:,2);
input3 = amplitude_of_peaks.after_ip_5_40_minutes (:,3);

figure
bar_plots_open_field({input1,input2,input3},{LegendVeh1;LegendDETQ;LegendVeh},{Color_Veh1;Color_DETQ;Color_Veh2}, name_plot,'Amplitude Ratio',fontSize_4figures, 1)
save_plots(name_plot)

%% 
save ('variables_5mice_average')