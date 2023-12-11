% Track mouse behavior from the Images 
% Author: Maria Chernysheva. Modified from "MouseActivity5" by Renzhi Han
% from Affiliation: The Ohio State University

clc; % clear the comman window
close all; % close all figures
clear; % erase all existing variables
MouseName='m033';
Folder='Y:\Chernysheva Maria\Fiber Photometry data\m033\20220723\OpenField_Veh2\182539320';
threshold=0.1;
MouseSize=100;

cd(Folder)
FolderFrames=[Folder,'\frames'];
% LogsFileName=[FolderLogs,'\logs'];
ListFrames=dir(FolderFrames); 
nframes=size(ListFrames,1)-2;

TimeStampFirstFrame_Date=ListFrames(4).name;
TimeStampLastFrame_Date=ListFrames(end).name;

Tstart= datevec(TimeStampFirstFrame_Date(length(TimeStampFirstFrame_Date)-20:length(TimeStampFirstFrame_Date)-4),'yyyymmddHHMMSSFFF');
Tend= datevec(TimeStampLastFrame_Date(length(TimeStampFirstFrame_Date)-20:length(TimeStampFirstFrame_Date)-4),'yyyymmddHHMMSSFFF');
SessionTime=etime(Tend,Tstart);
frameRate=nframes/SessionTime;

cd (FolderFrames)
I=imread([FolderFrames,'\',ListFrames(10000).name]);
[BW, xi, yi]=roipoly(I);
PositionROI=[xi(1) yi(1);xi(2) yi(2);xi(3) yi(3);xi(4) yi(4)];
xmin=min(PositionROI(:,1));
ymin=min(PositionROI(:,2));
xmax=max(PositionROI(:,1));
ymax=max(PositionROI(:,2));
width=xmax-xmin;
height=ymax-ymin;
cd(Folder)
%% for DETQ Open field project
jj=1;
Time_StartAnalyse=26; %*60; %in minutes, depending on when the OF started after the ip (can vary for each mouse
Time_EndAnalyse=Time_StartAnalyse+95; %*60;


Tstart_str=datestr(Tstart);
Tstart_num=datenum(Tstart);

Tstart_Analize_num=addtodate(Tstart_num,Time_StartAnalyse,'minute');
Tfinish_Analize_num=addtodate(Tstart_num,Time_EndAnalyse,'minute');


Tstart_Analize_vec=datevec(Tstart_Analize_num);
Tfinish_Analize_vec=datevec(Tfinish_Analize_num);

% 

for ii=1:nframes-1
Frames_Vec{ii}=datevec(ListFrames(ii+3).name(length(TimeStampFirstFrame_Date)-20:length(TimeStampFirstFrame_Date)-4),'yyyymmddHHMMSSFFF');
dT_startAnalyse(ii)=etime(Tstart_Analize_vec,Frames_Vec{ii});
dT_finishAnalyse(ii)=etime(Tfinish_Analize_vec,Frames_Vec{ii});
end

[m,Frame_Start]=min(abs(dT_startAnalyse));
[m,Frame_End]=min(abs(dT_finishAnalyse));

I=imread([FolderFrames,'\',ListFrames(Frame_Start).name]);
figure
imshow(I)
h = images.roi.Polygon(gca,'Position',PositionROI);
M = ~h.createMask();

% 
cd(Folder)
close all; % close all figures
% initiate the figure
  f = waitbar(0,'1','Name','Tracking mouse...');

MouseN=1;
m=1;
c1 = cell(MouseN, Frame_End-Frame_Start+1, 2);
     m1_majl = cell(MouseN, Frame_End-Frame_Start+1, 1);
     m1_minl = cell(MouseN, Frame_End-Frame_Start+1, 1);
     m1_ori = cell(MouseN, Frame_End-Frame_Start+1, 1);
     m1_ecc = cell(MouseN, Frame_End-Frame_Start+1, 1);
     m1_area = cell(MouseN, Frame_End-Frame_Start+1, 1);
     
for jj=Frame_Start:Frame_End
CurrentFrameName=ListFrames(jj).name; % take each frame 
   
CurrentFrameName2=[FolderFrames,'\',CurrentFrameName];
 
CurrentFrame=imread(CurrentFrameName2);

CurrentFrame(M) = 0;

A=rgb2gray(CurrentFrame);

 II=imbinarize(A,threshold);
 JJ=1-II;
 JJ(M) = 0;
 QQ=bwareaopen(JJ,MouseSize);
 
 QQ_crop = imcrop(QQ,[xmin ymin width height]);


 mouse= regionprops(QQ_crop,'Area','Centroid','MajorAxisLength','MinorAxisLength','Orientation','Eccentricity');
  if ~isempty([mouse.Area])                    
                    areaArray = [mouse.Area];
                    [~,idx] = max(areaArray);
                    c1{m}(jj,:) = floor(mouse(idx).Centroid +[xmin ymin]);
                    
                    m1_majl{m}(jj,:) = mouse(idx).MajorAxisLength;
                    m1_minl{m}(jj,:) = mouse(idx).MinorAxisLength;
                    m1_ori{m}(jj,:) = mouse(idx).Orientation;
                    m1_ecc{m}(jj,:) = mouse(idx).Eccentricity;
                    m1_area{m}(jj,:) = mouse(idx).Area;                    
                else
                    %the following is needed because mice sometimes jump
                    %up at the edge (out of the ROI).
                    c1{m}(jj,:) = c1{MouseN}(jj-1,:);
                    m1_majl{m}(jj,:) = m1_majl{m}(jj-1);
                    m1_minl{m}(jj,:) = m1_minl{m}(jj-1);
                    m1_ori{m}(jj,:) = m1_ori{m}(jj-1);
                    m1_ecc{m}(jj,:) = m1_ecc{m}(jj-1);
                    m1_area{m}(jj,:) = m1_area{m}(jj-1);
                 end
 
 

 waitbar((jj-Frame_Start)/(Frame_End-Frame_Start),f,sprintf(['Tracking Progress: ' num2str(round((jj-Frame_Start)*10000/(Frame_End-Frame_Start))/100) '%%']))
end

result.positions = {};
     result.area = {};
     result.orientation = {};
result.positions = [result.positions c1{m}];
         result.area = [result.area m1_area{m}];
         result.orientation = [result.orientation [m1_majl{m} m1_minl{m} m1_ori{m} m1_ecc{m}]];  
 save('result');



%% Read and pre-analise photometry recordings in Doric format 
%Note the .doric format is the same as .h5
%Author: Maria Chernysheva

%use cd function to chose a folder with data

% Folder='Y:\Chernysheva Maria\Fiber Photometry data\m022\20211103\OpenFieldVeh_Day3\144343427';
% give the name of the recordered file:
cd (Folder)

list_of_files=dir ('*.doric');
filename=list_of_files.name;

%get information about the Doric file:
DoricInfo = h5info(filename);

%read data recorded from defferent channels

%example line for reading a single channel:
%data = h5read(filename,'/Traces/Console/AIn-1 - Dem (AOut-1)/AIn-1 - Dem (AOut-1)' );

%create cells with datasets
%first colum: name of the channel
for ii=1:size (DoricInfo.Groups.Groups.Groups,1)
    DoricData{ii,1}=DoricInfo.Groups.Groups.Groups(ii).Datasets.Name ;
    CurrentChannel=sprintf('%s%s%s', DoricInfo.Groups.Groups.Groups(ii).Name,'/',DoricInfo.Groups.Groups.Groups(ii).Datasets.Name );
    DoricData{ii,2}=h5read(filename,CurrentChannel);
end
% 
% plot the raw demodulated data from two channels
%check that channel's names in the DoricData
%if the data were collected with the Program "UnpredictedReward":
%AI1:isosbestic - ch1
%AI2: functional -ch2
%DIO1-licking response - ch6
%DIO2-reward delivery - ch7
%DIO3-video frames - ch8
%DIO4-trigger doric- ch9
%Console time - ch9
FunctionalSignal_demodulated_channel=2;
IsosbesticSignal_demodulated_channel=1;
LickingResponse_channel=6;
VideoFrames_channel=8;
ConsoleTime_channel=9;

%extract the data to plot
FunctionalSignal_demodulated=DoricData{FunctionalSignal_demodulated_channel, 2};
IsosbesticSignal_demodulated=DoricData{IsosbesticSignal_demodulated_channel, 2};
LickingResponse=DoricData{LickingResponse_channel, 2};
VideoFrames_TTL=DoricData{VideoFrames_channel, 2};
ConsoleTime=DoricData{ConsoleTime_channel, 2};

%choose colors:
Color_FunctionalSignal=[0.4660 0.6740 0.1880];
Color_Licking=[0 0.4470 0.7410];

%Create a Plot Name based on the doric file
NamePlot=['RawSignals_', filename(1:end-6)];

figure
plot(ConsoleTime,FunctionalSignal_demodulated,'Color',Color_FunctionalSignal);
hold on
plot(ConsoleTime,IsosbesticSignal_demodulated,'k');
hold on

title(NamePlot,'Color','k','Fontsize',20, 'Interpreter', 'none')

ylabel('a.u.','Fontsize',26,'Color','k');
xlabel('time, s','Fontsize',26,'Color','k');
set(gca,'TickDir','out');
print(gcf, NamePlot,'-dpdf') 
savefig([NamePlot,'.fig'])
%% 

% clean the Doric signals from the Nan and big drops of the signal
ConsoleTime_Clean=ConsoleTime;
FunctionalSignal_demodulated_Clean=FunctionalSignal_demodulated;
IsosbesticSignal_demodulated_Clean=IsosbesticSignal_demodulated;
LickingResponse_Clean=LickingResponse;
VideoFrames_TTL_clean=VideoFrames_TTL;
DropSignalIndex2Delete=[];
k=1;
q=1;
for ii=1:length(IsosbesticSignal_demodulated)
if isnan(FunctionalSignal_demodulated(ii))==1||isnan(IsosbesticSignal_demodulated(ii))==1   % find Nan elements
    NaNIndex2Delete(k)=ii;
    k=k+1;
end

artifact_amplitude = 0.11;
if IsosbesticSignal_demodulated(ii)<artifact_amplitude||FunctionalSignal_demodulated(ii)<artifact_amplitude
    DropSignalIndex2Delete(q)=ii;
    q=q+1;
end
end

AllIndex2Delete=horzcat(NaNIndex2Delete,DropSignalIndex2Delete);

FunctionalSignal_demodulated_Clean(AllIndex2Delete)=[];
IsosbesticSignal_demodulated_Clean(AllIndex2Delete)=[];
ConsoleTime_Clean(AllIndex2Delete)=[];
LickingResponse_Clean(AllIndex2Delete)=[];
VideoFrames_TTL_clean(AllIndex2Delete)=[];

% Substract the isosbestic signal 

T_Baseline_Start=100;
T_BaselineFinish=1300;
dT_Doric_average=max(ConsoleTime)/length(FunctionalSignal_demodulated);

X=IsosbesticSignal_demodulated_Clean(round(T_Baseline_Start/dT_Doric_average):round(T_BaselineFinish/dT_Doric_average));
Y=FunctionalSignal_demodulated_Clean(round(T_Baseline_Start/dT_Doric_average):round(T_BaselineFinish/dT_Doric_average));
bls = polyfit(X,Y,1);

Y_fit_all = bls(1) .* IsosbesticSignal_demodulated_Clean + bls(2);
Y_dF_all = FunctionalSignal_demodulated_Clean - Y_fit_all;

dFF=100*((FunctionalSignal_demodulated_Clean - Y_fit_all)./Y_fit_all);
% 
cd(Folder)
figure

plot(ConsoleTime_Clean,dFF,'Color',Color_FunctionalSignal);
% NamePlot=['MotionCorrected V2 Base', filename(1:end-6)];
NamePlot=['MotionCorrected ', filename(1:end-6)];
title(NamePlot,'Color','k','Fontsize',20, 'Interpreter', 'none')

ylabel('dFF','Fontsize',26,'Color','k');
xlabel('time, s','Fontsize',26,'Color','k');
set(gca,'TickDir','out');
set(gca,'FontSize',16)
print(gcf, NamePlot,'-dpdf') 
savefig([NamePlot,'.fig'])
%


%% 
% Whole session Track
cd(Folder)
% MouseName='m033';
positions = result.positions;
orientation= result.orientation;
[x_center,y_center] = centroid(polyshape(xi,yi));
 polyout = scale(polyshape(xi',yi'),0.5,[x_center y_center]);
      MouseX=positions{1}(Frame_Start:Frame_End,1);
    MouseY=positions{1}(Frame_Start:Frame_End,2);
    
    MouseInCenter = inpolygon(MouseX,MouseY,polyout.Vertices(:,1),polyout.Vertices(:,2));
    timesMouseInCenter = numel(MouseX(MouseInCenter));
    timesMouseOutOfCenter = numel(MouseX(~MouseInCenter));
   
    
    FrameStep=SessionTime/nframes;
    TimeMouseInCenter=timesMouseInCenter*FrameStep; % in seconds
    TimeMouseOutOfCenter=timesMouseOutOfCenter*FrameStep; % in seconds
    
    Thigmotaxis_Percent=100*TimeMouseOutOfCenter/(TimeMouseOutOfCenter+TimeMouseInCenter);
    TimePercent_CenterZone=100-Thigmotaxis_Percent;


figure

plot(positions{1}(Frame_Start:Frame_End,1),positions{1}(Frame_Start:Frame_End,2),'Color','k');
hold on
% rectangle('Position',[xmin,ymin, width, height])
plot(xi,yi,'k','LineWidth',2); % plot open field arena for each mouse
hold on
plot(polyout) %,'--','k','LineWidth',2);
hold on

plot(MouseX(MouseInCenter),MouseY(MouseInCenter),'r.','MarkerSize',5)
set ( gca, 'ydir', 'reverse' )
axis ('square')
axis off
NamePlot=['Mouse Track ', MouseName]; %,filename(1:end-6)];

% title(['Whole Session Track. ', 'Thigmotaxis=',num2str(Thigmotaxis_Percent),'%'; {filename(1:end-6)}],'Color','k','Fontsize',20, 'Interpreter', 'none')
% set (gcf, 'Position',[9.15458333333333                         5          15.8454166666667                 16.034375])
title(['Time in the center=',num2str(TimePercent_CenterZone),'%'; {MouseName}],'Color','k','Fontsize',20, 'Interpreter', 'none')
print(gcf, NamePlot,'-dpdf') 
savefig([NamePlot,'.fig'])

positions = result.positions;
orientation= result.orientation;
area = result.area;
   dd1 = 0;
   distance1 = [0];
   pathl1 =[0];
mean_pixel_scale=1/(width+height); %in meters

nframes_Video=Frame_End-Frame_Start+1;

TimeStampFirstFrame_Video=ListFrames(Frame_Start).name;
TimeStampLastFrame_Video=ListFrames(Frame_End).name;

Tstart_Video= datevec(TimeStampFirstFrame_Video(length(TimeStampFirstFrame_Video)-20:length(TimeStampFirstFrame_Video)-4),'yyyymmddHHMMSSFFF');
Tend_Video= datevec(TimeStampLastFrame_Video(length(TimeStampFirstFrame_Video)-20:length(TimeStampFirstFrame_Video)-4),'yyyymmddHHMMSSFFF');
VideoTime=etime(Tend_Video,Tstart_Video);

FrameStep=VideoTime/nframes_Video;

        for k = Frame_Start+1:size(positions{1},1)
              D1 = sqrt((positions{1}(k,1)-positions{1}(k-1,1))^2+(positions{1}(k,2)-positions{1}(k-1,2))^2);
              D1 = D1.* mean_pixel_scale; % path length in mm
              pathl1 = [pathl1; D1];
              dd1 = dd1 + D1; % travel distance in m
              distance1 = [distance1; dd1];
        end
        
        time_distance1=0:FrameStep:VideoTime-FrameStep;
        
        figure
        plot(time_distance1,distance1,'Color','k','LineWidth',2)
        xlabel('Time (s)');
    ylabel('Traveled Distance (m)');
    
% NamePlot=['Locomotion ',filename(1:end-6)]; for NAc open field
NamePlot=['Locomotion ',MouseName];

% title(['Whole Session Locomotion. ', 'Distance=',num2str(max(distance1)),'m'; {filename(1:end-6)}],'Color','k','Fontsize',20, 'Interpreter', 'none')

title(['Distance=',num2str(max(distance1)),'m'; {MouseName}],'Color','k','Fontsize',20, 'Interpreter', 'none') %for PD mice

print(gcf, NamePlot,'-dpdf') 
savefig([NamePlot,'.fig'])

% average Thigmotaxis and Distance per 10 minutes (or other loop)

kk=5;
jj=1;
close all
for kk=5:5:40
    
    LoopMinutes=kk;
    
LoopSeconds=LoopMinutes*60;
NLoops=100/LoopMinutes; %2 hours recording , 120 minutes
k=2;
FrameEndLoop=[];
distance_Loop=[];

nFramesLoop=round(LoopSeconds/FrameStep);


%Option 2: find exactly the Frame Loop of 10 minutes by the timesteps
T_startLoop=[];
T_startLoop(1)=Time_StartAnalyse;
Tstart_num=datenum(Tstart);
Frame_StartLoop=[];
for k=2:NLoops
    T_startLoop(k)=T_startLoop(k-1)+LoopMinutes;
    T_startLoop_num(k)=addtodate(Tstart_num,T_startLoop(k),'minute');
    T_startLoop_vec{k}=datevec(T_startLoop_num(k));

dT_startLoop=[];

for ii=1:nframes-1
% Frames_Vec{ii}=datevec(ListFrames(ii+2).name(length(TimeStampFirstFrame_Date)-20:length(TimeStampFirstFrame_Date)-4),'yyyymmddHHMMSSFFF');
dT_startLoop(ii)=etime(T_startLoop_vec{k},Frames_Vec{ii});

end

[m,Frame_StartLoop(k)]=min(abs(dT_startLoop));

end


Frame_StartLoop(1)=Frame_Start;
Frame_StartLoop=[Frame_StartLoop,Frame_End];
Thigmotaxis_Percent_Loop=[];

k=1;
for k=1:NLoops
    MouseX=positions{1}(Frame_StartLoop(k):Frame_StartLoop(k+1),1);
    MouseY=positions{1}(Frame_StartLoop(k):Frame_StartLoop(k+1),2);
    
    MouseInCenter = inpolygon(MouseX,MouseY,polyout.Vertices(:,1),polyout.Vertices(:,2));
    timesMouseInCenter = numel(MouseX(MouseInCenter));
    timesMouseOutOfCenter = numel(MouseX(~MouseInCenter));
   
    

    TimeMouseInCenter=timesMouseInCenter*FrameStep; % in seconds
    TimeMouseOutOfCenter=timesMouseOutOfCenter*FrameStep; % in seconds
    
    Thigmotaxis_Percent_Loop(k)=100*TimeMouseOutOfCenter/(TimeMouseOutOfCenter+TimeMouseInCenter);
    
    distance_Loop(k)=distance1(Frame_StartLoop(k+1)-Frame_Start+1)-distance1(Frame_StartLoop(k)-Frame_Start+1);
end

Thigmotaxis_perLoop{2,jj}=Thigmotaxis_Percent_Loop;
distance_perLoop{2,jj}=distance_Loop;

Thigmotaxis_perLoop{1,jj}=[num2str(kk),' min'];
distance_perLoop{1,jj}=[num2str(kk),' min'];

timeLoop=[];
timeLoop=LoopMinutes:LoopMinutes:NLoops*LoopMinutes;

figure
plot(timeLoop,Thigmotaxis_Percent_Loop,'.','Color','k','MarkerSize',10)
ylabel('Average Thigmotaxis (%)');
xlabel('Time, min');
axis([-10 130 0 100])
NamePlot=['Average Thigmotaxis ',num2str(LoopMinutes),' min ',filename(1:end-6)];

title(['Average Thigmotaxis. ', 'Thigmotaxis per ',num2str(LoopMinutes),' min'; {filename(1:end-6)}],'Color','k','Fontsize',20, 'Interpreter', 'none')
print(gcf, NamePlot,'-dpdf') 
savefig([NamePlot,'.fig'])
pause(0.4)

figure
plot(timeLoop,distance_Loop,'.','Color','k','MarkerSize',10)
ylabel('Traveled Distance (m)');
xlabel('Time, min');
axis([-10 130 0 140])
NamePlot=['Average Locomotion ',num2str(LoopMinutes),' min ',filename(1:end-6)];

title(['Average Locomotion. ', 'Distance per ',num2str(LoopMinutes),' min'; {filename(1:end-6)}],'Color','k','Fontsize',20, 'Interpreter', 'none')
print(gcf, NamePlot,'-dpdf') 
savefig([NamePlot,'.fig'])

jj=jj+1;
end


%% 

% save all variables

close all
save ('AllVariables')
% 
%% 
cd(Folder)
varList=who;
StructVar_name='m034_OpenField_DETQ_Day2.mat';
StructVar=struct;
%use dynamic fieldnames

for index = 1:numel(varList)
  
    StructVar.(varList{index}) = eval(varList{index});
end 

m034_OpenField_DETQ_Day2=StructVar;
save (StructVar_name,'m034_OpenField_DETQ_Day2')

