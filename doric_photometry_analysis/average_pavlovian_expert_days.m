%% averaging expert days

%load variables

close all
clearvars
clc

%give folder with matlab functions
folder_with_matlab_scripts = 'Y:\Chernysheva Maria\MATLAB\Matlab_scripts_AlloLight';
addpath (genpath(folder_with_matlab_scripts))

data_folder.Veh1 = 'Y:\Chernysheva Maria\Fiber_Photometry_data\data_for_averaging\Pavlovian\PreTrialFit_Baseline\data_to_average\veh_day1';
data_folder.DETQ = 'Y:\Chernysheva Maria\Fiber_Photometry_data\data_for_averaging\Pavlovian\PreTrialFit_Baseline\data_to_average\detq';
data_folder.Veh2 = 'Y:\Chernysheva Maria\Fiber_Photometry_data\data_for_averaging\Pavlovian\PreTrialFit_Baseline\data_to_average\veh_day2';
Folder_4Figures = 'Y:\Chernysheva Maria\Fiber_Photometry_data\data_for_averaging\Pavlovian\PreTrialFit_Baseline\data_to_average\figures';
%%
number_of_different_groups = 3;
folders = fieldnames(data_folder);
RangeTrial=19;
for group_number = 1:number_of_different_groups
    
    cd  (data_folder.(folders{group_number}))
    list_of_files = dir ('*.mat'); %list the data folder content
    
    for mouse_Number = 1:1:size(list_of_files,1) %loop mice
        name_file{mouse_Number} = list_of_files(mouse_Number).name; %Find the names of the files for each mouse
        load(name_file{mouse_Number}) %load the file
        data_injection_type{group_number}{mouse_Number} = StrStr;
        name_current_file = fieldnames(data_injection_type{group_number}{mouse_Number});
        current_deltaFF.(folders{group_number}){mouse_Number} = data_injection_type{1, group_number}{mouse_Number}.(char(name_current_file)).deltaFF_RewardedSound(1:RangeTrial,:);
        clearvars StrStr
    end
    deltaFF.(folders{group_number}) = vertcat (current_deltaFF.(folders{group_number}){:});
end

%% average dFF signal
SecondsBefore = 5;
SecondsAfter = 15;
movWindow = 20;
NMice = 6;


LineWidth = 1.2;

Color_Grey = [0.7 0.7 0.7];
Color_NiceRed = [255, 128, 128]./255;
Color_Veh1 = [0.2 0.2 0.2];
Color_DETQ = Color_NiceRed;
Color_Veh2 = Color_Grey;
Color.(folders{1}) = Color_Veh1;
Color.(folders{2}) = Color_DETQ;
Color.(folders{3}) = Color_Veh2;

LegendVeh1=['\color[rgb]{',num2str(Color_Veh1),'}Veh'];
LegendDETQ=['\color[rgb]{',num2str(Color_DETQ),'}DETQ'];
LegendVeh=['\color[rgb]{',num2str(Color_Veh2),'}Veh'];

dT_Doric_average=(SecondsBefore+SecondsAfter)/(size(deltaFF.(folders{group_number}),2)-1);

step_before=round(SecondsBefore/dT_Doric_average);
index_4_seconds = round(4/dT_Doric_average);
index_delay = round(10.4/dT_Doric_average);
step_after=round(SecondsAfter/dT_Doric_average);
Time2Reward=-1*(SecondsBefore):dT_Doric_average:SecondsAfter;



for group_number = 1:number_of_different_groups
    
    for each_trial = 1:size (deltaFF.(folders{group_number}),1)
        movmean_deltaFF.(folders{group_number})(each_trial,:) = movmean(deltaFF.(folders{group_number})(each_trial,:),movWindow);
    end
    
    Mean_deltaFF.(folders{group_number}) = mean(movmean_deltaFF.(folders{group_number}),1);
    sem_Mean_deltaFF.(folders{group_number}) = std(movmean_deltaFF.(folders{group_number}),1)/sqrt(size(movmean_deltaFF.(folders{group_number}),1));
    
end

deltaT_CS_US =10.4551;

%% average trace, Figure 7e

figure
NamePlot='All mice Expert Days MovMean Design2';
tickLength1=0.025;
tickLength2=0.1;

for group_number = 1:number_of_different_groups
H.(folders{group_number}) = shadedErrorBar(Time2Reward,Mean_deltaFF.(folders{group_number}),sem_Mean_deltaFF.(folders{group_number}),{'color', Color.(folders{group_number}),'linewidth',LineWidth},0);
main_lines{group_number} = H.(folders{group_number}).mainLine;
hold on
end

xline(deltaT_CS_US,':')

title ([NamePlot,{''}])
hold on
xline(0,':')
hold on
axis([min(Time2Reward) max(Time2Reward) -1 9])

ylabel('\DeltaF/F (%)','Fontsize',14,'Color','k');
xlabel('Time from CS onset, s','Fontsize',14,'Color','k');
title('Expert Sessions')

legend(vertcat(main_lines{:}),{LegendVeh1, ['\color[rgb]{',num2str(Color_DETQ),'}DETQ'], ['\color[rgb]{',num2str(Color_Veh2),'}Veh']},'Location','eastoutside','FontSize',16)

set(gca,'FontName','Arial');
set(gca,'FontSize',16);
set(gca,'TickDir','out');
set(gca,'TickLength',[tickLength1, tickLength2])

box off
legend boxoff

savefig([NamePlot,'.fig'])

%% Fig7f left:
%Plot design:

x1=0.8;
x2=1.8;
x3=2.8;
GreyLineWidth=0.5;
DataMarkerSize=8;
ErrorBarWidth=12;
ErrorBarLineWidth=1.5;
% proportion of the figure:
pbaspect_1=1;
pbaspect_2=2;
pbaspect_3=1;
step_before=round(SecondsBefore/dT_Doric_average);

for group_number = 1:number_of_different_groups
    
    for mouse_Number = 1:1:size(list_of_files,1)
name_current_file = fieldnames(data_injection_type{group_number}{mouse_Number});
current_peak_dFF.(folders{group_number}){mouse_Number} = max(mean(data_injection_type{1, group_number}{mouse_Number}.(char(name_current_file)).deltaFF_RewardedSound(1:RangeTrial,step_before:step_before+index_4_seconds),1));
    end
    
    peak_dFF_per_group.(folders{group_number}) = vertcat(current_peak_dFF.(folders{group_number}){:});
    
end

name_plot='Peak CS CSBaseline';
%
LegendVeh1=['\color[rgb]{',num2str(Color_Veh1),'}Veh'];
LegendDETQ=['\color[rgb]{',num2str(Color_DETQ),'}DETQ'];
LegendVeh=['\color[rgb]{',num2str(Color_Veh2),'}Veh'];


input1 = peak_dFF_per_group.(folders{1});
input2 = peak_dFF_per_group.(folders{2});
input3 = peak_dFF_per_group.(folders{3});

fontSize_4figures = 12;

figure
SEM_Input = bar_plots_open_field({input1,input2,input3},{LegendVeh1;LegendDETQ;LegendVeh},{Color_Veh1;Color_DETQ;Color_Veh2}, name_plot,'Peak\DeltaF/F_0 (%)',fontSize_4figures, 1);

save_plots(name_plot)

excel_source_data_fig7_panel_f_CS = vertcat ([input1, input2, input3],...
    [mean(input1),mean(input2),mean(input3)],...
    SEM_Input);

%% plot the heat maps. Fig 7d
mouse_number = 2;
mouse_name = 'm011';
NamePlot=['Heat Map 8 Design2 ', mouse_name];

number_of_trials = 19;
first_trial_per_mouse = (mouse_number-1)*number_of_trials+1;

clims = [0 10];
tickLength1 = 0.025;
tickLength2 = 0.1;
F = figure;

range_2plot=19+19; 

C=movmean_deltaFF.(folders{1});
x=Time2Reward;
y=(1:size(movmean_deltaFF.(folders{1})(first_trial_per_mouse:first_trial_per_mouse+number_of_trials-1,:),1));

subplot(1,3,1)

imagesc(x,y,movmean_deltaFF.(folders{1})(first_trial_per_mouse:first_trial_per_mouse+number_of_trials-1,:),clims);
xlabel('Time from CS onset, s', 'FontSize', 16);
hold on
xline(0,'w','Linewidth',2)
hold on
xline(deltaT_CS_US,'w','Linewidth',2)

title(['Veh',{''}], 'FontSize', 16,'Color',Color_Veh1);
set(gca,'FontName','Arial');
set(gca,'FontSize',16);
ylabel('Trials', 'FontSize', 20);
set(gca,'TickDir','out');
set(gca,'TickLength',[tickLength1, tickLength2])
box off

hold on
s2=subplot(1,3,2);
imagesc(x,y,movmean_deltaFF.(folders{2})(first_trial_per_mouse:first_trial_per_mouse+number_of_trials-1,:),clims);
box off
hold on
xline(0,'w','Linewidth',2)
hold on
xline(deltaT_CS_US,'w','Linewidth',2)
title(['DETQ',{''}],'Color',Color_DETQ);
xlabel('Time from CS onset, s', 'FontSize', 16);
set(gca,'FontName','Arial');
% set(gca,'FontSize',16);
set(gca,'TickDir','out');
set(gca,'TickLength',[tickLength1, tickLength2])
box off
set(gca,'Yticklabel',[])

hold on
s3=subplot(1,3,3);
imagesc(x,y,movmean_deltaFF.(folders{3})(first_trial_per_mouse:first_trial_per_mouse+number_of_trials-1,:),clims);
hold on
xline(0,'w','Linewidth',2)
hold on
xline(deltaT_CS_US,'w','Linewidth',2)
title(['Veh',{''}],'Color',Color_Veh2);
xlabel('Time from CS onset, s', 'FontSize', 16);
box off
set(gca,'Yticklabel',[])

set(gca,'FontName','Arial');
set(gca,'FontSize',16);
set(gca,'TickDir','out');
set(gca,'TickLength',[tickLength1, tickLength2])

h=colorbar;
h_ylablel=ylabel(h,'\DeltaF/F (%)','Rotation',270.0,'FontSize',20);
h.Label.Position(1) = 4;
set(gcf, 'Position', get(0, 'Screensize'))
set(gcf,'PaperOrientation','landscape');

set(s2,'Position',[0.3808    0.1100    0.2134    0.8150],'FontSize',16);
set(s3,'Position',[0.6316    0.1100    0.2134    0.8150],'FontSize',16);

%% Fig 7f right
%Reward

for group_number = 1:number_of_different_groups
    
    for mouse_Number = 1:1:size(list_of_files,1)
name_current_file = fieldnames(data_injection_type{group_number}{mouse_Number});
current_peak_dFF_UC.(folders{group_number}){mouse_Number} = max(mean(data_injection_type{1, group_number}{mouse_Number}.(char(name_current_file)).deltaFF_Input_Reward(1:RangeTrial,step_before:step_before+index_4_seconds),1));
    end
    
    peak_dFF_UC_per_group.(folders{group_number}) = vertcat(current_peak_dFF_UC.(folders{group_number}){:});
    
end




name_plot='Peak US';
LegendVeh1=['\color[rgb]{',num2str(Color_Veh1),'}Veh'];
LegendDETQ=['\color[rgb]{',num2str(Color_DETQ),'}DETQ'];
LegendVeh=['\color[rgb]{',num2str(Color_Veh2),'}Veh'];

input1 = peak_dFF_UC_per_group.(folders{1});
input2 = peak_dFF_UC_per_group.(folders{2});
input3 = peak_dFF_UC_per_group.(folders{3});
figure
SEM_Input = bar_plots_open_field({input1,input2,input3},{LegendVeh1;LegendDETQ;LegendVeh},{Color_Veh1;Color_DETQ;Color_Veh2}, name_plot,'Peak\DeltaF/F_0 (%)',fontSize_4figures, 1);

save_plots(name_plot)

%% Latensy to lick
%!don't forget to clean the latensy from the mistakes (e.g. latensy more
%then 12 seconds (the duration of the sound)
data_injection_type{1, 1}{1}.(char(fieldnames(data_injection_type{1}{1}))).Latency_Sound_Reward_Time(14) = [];

for group_number = 1:number_of_different_groups
    
    for mouse_Number = 1:1:size(list_of_files,1)
name_current_file = fieldnames(data_injection_type{group_number}{mouse_Number});
current_Latency_Sound_Reward_Time.(folders{group_number}){mouse_Number} = mean(data_injection_type{1, group_number}{mouse_Number}.(char(name_current_file)).Latency_Sound_Reward_Time(1:RangeTrial));
    end
    
    Latency_Sound_Reward_Time_per_group.(folders{group_number}) = vertcat(current_Latency_Sound_Reward_Time.(folders{group_number}){:});
    
end

name_plot='Latency to Lick';
LegendVeh1=['\color[rgb]{',num2str(Color_Veh1),'}Veh'];
LegendDETQ=['\color[rgb]{',num2str(Color_DETQ),'}DETQ'];
LegendVeh=['\color[rgb]{',num2str(Color_Veh2),'}Veh'];

input1 = Latency_Sound_Reward_Time_per_group.(folders{1});
input2 = Latency_Sound_Reward_Time_per_group.(folders{2});
input3 = Latency_Sound_Reward_Time_per_group.(folders{3});
figure
SEM_Input = bar_plots_open_field({input1,input2,input3},{LegendVeh1;LegendDETQ;LegendVeh},{Color_Veh1;Color_DETQ;Color_Veh2}, name_plot,'Latency to Lick',fontSize_4figures, 1);

save_plots(name_plot)

%% Proportion of licking during sound Fig7g right
for group_number = 1:number_of_different_groups
    
    for mouse_Number = 1:1:size(list_of_files,1)
name_current_file = fieldnames(data_injection_type{group_number}{mouse_Number});
current_LickPercentProportion.(folders{group_number}){mouse_Number} = lick_proportion(data_injection_type{1, group_number}{mouse_Number}.(char(name_current_file)).LickVector_sound(1:RangeTrial,:),...
                                                                                      data_injection_type{1, group_number}{mouse_Number}.(char(name_current_file)).LickVector_NoSound);
    end
    
    LickPercentProportion_per_group.(folders{group_number}) = vertcat(current_LickPercentProportion.(folders{group_number}){:});
    
end

name_plot = 'Licking proportion';
LegendVeh1=['\color[rgb]{',num2str(Color_Veh1),'}Veh'];
LegendDETQ=['\color[rgb]{',num2str(Color_DETQ),'}DETQ'];
LegendVeh=['\color[rgb]{',num2str(Color_Veh2),'}Veh'];

input1 = LickPercentProportion_per_group.(folders{1});
input2 = LickPercentProportion_per_group.(folders{2});
input3 = LickPercentProportion_per_group.(folders{3});
figure
SEM_Input = bar_plots_open_field({input1,input2,input3},{LegendVeh1;LegendDETQ;LegendVeh},{Color_Veh1;Color_DETQ;Color_Veh2}, name_plot,['Licking during', {'CS vs. No CS'}],fontSize_4figures, 1);

save_plots(name_plot)

% excel_source_data_Fig6_panelg_right  = vertcat ([input1, input2, input3],...
%     [mean(input1),mean(input2),mean(input3)],...
%     [SEM_Input(1), SEM_Input(2), SEM_Input(3)]);
