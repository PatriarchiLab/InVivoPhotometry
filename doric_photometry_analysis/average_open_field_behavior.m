%averaging Open field data
%Author: Maria Chernysheva

%load the data
close all
clearvars
clc
%give folder with matlab functions
folder_with_matlab_scripts = 'Y:\Chernysheva Maria\MATLAB\Matlab_scripts_AlloLight';
addpath (genpath(folder_with_matlab_scripts))

%give the folder with the data subfolders
main_folder_with_data_subfolders = 'Y:\Chernysheva Maria\Fiber_Photometry_data\data_for_averaging\Open_Field\dLight13b';

cd (main_folder_with_data_subfolders)
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 5 20 10]);% define the default size for the matlab figure

group_name{1} = 'veh_day1'; 
group_name{2} = 'detq_day2'; 
group_name{3} = 'veh_day3'; 

%indicate minutes from the start of the recording (recorded in the experimental log
%book)
time_ip_veh_day1 = [27,  24,  26,   24.5, 28]; %m008, m012, m020, m021, m022
time_ip_detq_day2 = [32, 34.5,   26.5, 24,   26]; %m008, m012, m020, m021, m022
time_ip_veh_day3 = [28.5, 27.5, 26.5, 31.5, 25]; %m008, m012, m020, m021, m022


%!uncomment for the AlloLite-ctr mice:
% time_ip_veh_day1=[33,25,27,25]; %m027, m028, m033, m034
% time_ip_detq_day2=[28,27,26,26]; %m027, m028, m033, m034
% time_ip_veh_day3=[29,28,26,25]; %m027, m028, m033, m034


index_ip.veh_day1 = round(time_ip_veh_day1*60/0.0082);
index_ip.detq_day2 = round(time_ip_detq_day2*60/0.0082);
index_ip.veh_day3 = round(time_ip_veh_day3*60/0.0082);


index_20minutes = round(20*60/0.0082);
index_90minutes = round(90*60/0.0082);
% 
number_of_different_groups = 3; 

for group_number=1:number_of_different_groups
folders_with_data{group_number} = uigetdir(main_folder_with_data_subfolders,['Select Folder with the Drug ',num2str(group_number)]);
end

%set the folder where you would like to store Figures and Variables for
%this experiment: 

folder_4_figures = uigetdir(main_folder_with_data_subfolders,'Select Folder to save Plots and Vars');

%Automatically create sub-folders to have dFF analysed classical way, dFF_465, dFF_405
cd (folder_4_figures)
%% 
for group_number = 1:number_of_different_groups %loop drugs (veh vs DETQ etc)
    
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

% 
fontSize_4figures=12;
Color_Grey=[0.7 0.7 0.7];
Color_NiceRed=[255, 128, 128]./255;

Color_Veh1=[0.2 0.2 0.2];
Color_DETQ=Color_NiceRed;
Color_Veh2=Color_Grey;

LegendVeh1=['\color[rgb]{',num2str(Color_Veh1),'}Veh'];
LegendDETQ=['\color[rgb]{',num2str(Color_DETQ),'}DETQ'];
LegendVeh=['\color[rgb]{',num2str(Color_Veh2),'}Veh'];

%% extended Figure 5f, left. Distance per 10 min after ip
interval = 2; %2: 10 minutes intervals
NIntervals = 9; % 9 intervals of 10 minutes within 90 minutes

for group_number = 1:number_of_different_groups
    
    for mouse_Number = 1:1:size(list_of_files,1)
name_current_file = fieldnames(data_cohort{group_number}{mouse_Number});
current_distance_per_loop.(group_name{group_number}){mouse_Number} = data_cohort{1, group_number}{mouse_Number}.(char(name_current_file)).distance_perLoop{2,interval}(1:NIntervals);
    end
    
    distance_per_loop_group.(group_name{group_number}) = vertcat(current_distance_per_loop.(group_name{group_number}){:});
    
end

name_plot='Distance per 10 min after ip';

input1 = distance_per_loop_group.(group_name{1});
input2 = distance_per_loop_group.(group_name{2});
input3 = distance_per_loop_group.(group_name{3});

SEM_Input_Day1=std(input1(:,1))/sqrt(length(input1(:,1)));
SEM_Input_Day2=std(input2(:,1))/sqrt(length(input2(:,1)));
SEM_Input_Day3=std(input3(:,1))/sqrt(length(input3(:,1)));

Mean_Day1=mean(input1(:,1));
Mean_Day2=mean(input2(:,1));
Mean_Day3=mean(input3(:,1));

for ii=2:size (input1,2)
SEM_Input_Day1=[SEM_Input_Day1 std(input1(:,ii))/sqrt(length(input1(:,ii)))];
SEM_Input_Day2=[SEM_Input_Day2 std(input2(:,ii))/sqrt(length(input2(:,ii)))];
SEM_Input_Day3=[SEM_Input_Day3 std(input3(:,ii))/sqrt(length(input3(:,ii)))];

Mean_Day1=[Mean_Day1 mean(input1(:,ii))];
Mean_Day2=[Mean_Day2 mean(input2(:,ii))];
Mean_Day3=[Mean_Day3 mean(input3(:,ii))];
end

figure
e1=errorbar(1:1:size (input1,2),Mean_Day1,SEM_Input_Day1(1:end),'o','LineWidth',1,'MarkerSize',8,'MarkerFaceColor',Color_Veh1);
e1.Color = Color_Veh1;
hold on
e2=errorbar(1:1:size (input1,2),Mean_Day2,SEM_Input_Day2(1:end),'o','LineWidth',1,'MarkerSize',8,'MarkerFaceColor',Color_DETQ);
hold on
e2.Color = Color_DETQ;
e3=errorbar(1:1:size (input1,2),Mean_Day3,SEM_Input_Day3(1:end),'o','LineWidth',1,'MarkerSize',8,'MarkerFaceColor',Color_Veh2);
e3.Color = Color_Veh2;
hold on
axis([0 (NIntervals+1) 15 35])
set(gca,'XTick',[0:1:(NIntervals+1)]);

title ([name_plot,{''}])
ylabel('Distance per 10 min, m','Fontsize',fontSize_4figures,'Color','k');
xlabel('10 min intervals','Fontsize',fontSize_4figures,'Color','k');

set(gca,'TickDir','out');
set(gca,'TickLength',[0.025, 0.1])
box off
savefig([name_plot,'.fig'])

%% Average Distance per 90 min. Extended Figure 5f, right
name_plot='Total distance after ip';

average_distance_day1 = sum(distance_per_loop_group.(group_name{1}),2);
average_distance_day2 = sum(distance_per_loop_group.(group_name{2}),2);
average_distance_day3 = sum(distance_per_loop_group.(group_name{3}),2);

input1 = average_distance_day1;
input2 = average_distance_day2;
input3 = average_distance_day3;

figure
SEM_Input = bar_plots_open_field({input1,input2,input3},{LegendVeh1;LegendDETQ;LegendVeh},{Color_Veh1;Color_DETQ;Color_Veh2}, name_plot,'Total Distance',fontSize_4figures, 1);
 
%% extended Figure 5g, left. Time in the center per 10 min after ip
interval = 2; %2: 10 minutes intervals
NIntervals = 9; % 9 intervals of 10 minutes within 90 minutes

for group_number = 1:number_of_different_groups
    
    for mouse_Number = 1:1:size(list_of_files,1)
name_current_file = fieldnames(data_cohort{group_number}{mouse_Number});
current_thigmotaxis_per_loop.(group_name{group_number}){mouse_Number} = data_cohort{1, group_number}{mouse_Number}.(char(name_current_file)).Thigmotaxis_perLoop{2,2}(1:9);
    end
    
    thigmotaxis_per_loop_group.(group_name{group_number}) = vertcat(current_thigmotaxis_per_loop.(group_name{group_number}){:});
    
end

name_plot='Time in the center per 10 min after ip';

input1 = 100-thigmotaxis_per_loop_group.(group_name{1});
input2 = 100-thigmotaxis_per_loop_group.(group_name{2});
input3 = 100-thigmotaxis_per_loop_group.(group_name{3});

SEM_Input_Day1=std(input1(:,1))/sqrt(length(input1(:,1)));
SEM_Input_Day2=std(input2(:,1))/sqrt(length(input2(:,1)));
SEM_Input_Day3=std(input3(:,1))/sqrt(length(input3(:,1)));

Mean_Day1=mean(input1(:,1));
Mean_Day2=mean(input2(:,1));
Mean_Day3=mean(input3(:,1));

for ii=2:size (input1,2)
SEM_Input_Day1=[SEM_Input_Day1 std(input1(:,ii))/sqrt(length(input1(:,ii)))];
SEM_Input_Day2=[SEM_Input_Day2 std(input2(:,ii))/sqrt(length(input2(:,ii)))];
SEM_Input_Day3=[SEM_Input_Day3 std(input3(:,ii))/sqrt(length(input3(:,ii)))];

Mean_Day1=[Mean_Day1 mean(input1(:,ii))];
Mean_Day2=[Mean_Day2 mean(input2(:,ii))];
Mean_Day3=[Mean_Day3 mean(input3(:,ii))];
end

figure
e1=errorbar(1:1:size (input1,2),Mean_Day1,SEM_Input_Day1(1:end),'o','LineWidth',1,'MarkerSize',8,'MarkerFaceColor',Color_Veh1);
e1.Color = Color_Veh1;
hold on
e2=errorbar(1:1:size (input1,2),Mean_Day2,SEM_Input_Day2(1:end),'o','LineWidth',1,'MarkerSize',8,'MarkerFaceColor',Color_DETQ);
hold on
e2.Color = Color_DETQ;
e3=errorbar(1:1:size (input1,2),Mean_Day3,SEM_Input_Day3(1:end),'o','LineWidth',1,'MarkerSize',8,'MarkerFaceColor',Color_Veh2);
e3.Color = Color_Veh2;
hold on
axis([0 (NIntervals+1) 0 30])
set(gca,'XTick',[0:1:(NIntervals+1)]);

title ([name_plot,{''}])
ylabel('Center time per 10 min, m','Fontsize',fontSize_4figures,'Color','k');
xlabel('10 min intervals','Fontsize',fontSize_4figures,'Color','k');

set(gca,'TickDir','out');
set(gca,'TickLength',[0.025, 0.1])
box off
savefig([name_plot,'.fig'])

%% Average time in the center per 90 minutes

name_plot='Total time in the Center after ip';

for group_number = 1:number_of_different_groups
    
    for mouse_Number = 1:1:size(list_of_files,1)
name_current_file = fieldnames(data_cohort{group_number}{mouse_Number});
current_center_time.(group_name{group_number}){mouse_Number} = time_in_the_center_90min_OF(data_cohort{1, group_number}{mouse_Number}.(char(name_current_file)));
    end
    
    TimePercent_CenterZone_group.(group_name{group_number}) = vertcat(current_center_time.(group_name{group_number}){:});
    
end

input1=TimePercent_CenterZone_group.(group_name{1});
input2=TimePercent_CenterZone_group.(group_name{2});
input3=TimePercent_CenterZone_group.(group_name{3});

figure
bar_plots_open_field({input1,input2,input3},{LegendVeh1;LegendDETQ;LegendVeh},{Color_Veh1;Color_DETQ;Color_Veh2}, name_plot,'Time in the Center, %',fontSize_4figures, 1)
SavePlots(name_plot)                                      
