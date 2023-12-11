%% Average time in the center per 90 minutes
function TimePercent_CenterZone_90min = time_in_the_center_90min_OF(structure)
% Load:

Time_StartAnalyse=structure.Time_StartAnalyse ;
Frame_Start=structure.Frame_Start;
Frames_Vec=structure.Frames_Vec;
xi=structure.xi;
yi=structure.yi;
result=structure.result;
TimeStampFirstFrame_Date=structure.TimeStampFirstFrame_Date;
nframes=size(structure.ListFrames,1)-3;
FrameStep=structure.FrameStep;

Tstart= datevec(TimeStampFirstFrame_Date(length(TimeStampFirstFrame_Date)-20:length(TimeStampFirstFrame_Date)-4),'yyyymmddHHMMSSFFF');
Time_EndAnalyse=Time_StartAnalyse+90; %*60;


Tstart_str=datestr(Tstart);
Tstart_num=datenum(Tstart);

Tstart_Analize_num=addtodate(Tstart_num,Time_StartAnalyse,'minute');
Tfinish_Analize_num=addtodate(Tstart_num,Time_EndAnalyse,'minute');


Tstart_Analize_vec=datevec(Tstart_Analize_num);
Tfinish_Analize_vec=datevec(Tfinish_Analize_num);

for ii=1:nframes
dT_finishAnalyse(ii)=etime(Tfinish_Analize_vec,Frames_Vec{ii});
end

[m,Frame_End]=min(abs(dT_finishAnalyse));

positions = result.positions;

[x_center,y_center] = centroid(polyshape(xi,yi));
 polyout = scale(polyshape(xi',yi'),0.5,[x_center y_center]);
      MouseX=positions{1}(Frame_Start:Frame_End,1);
    MouseY=positions{1}(Frame_Start:Frame_End,2);
    
    MouseInCenter = inpolygon(MouseX,MouseY,polyout.Vertices(:,1),polyout.Vertices(:,2));
    timesMouseInCenter = numel(MouseX(MouseInCenter));
    timesMouseOutOfCenter = numel(MouseX(~MouseInCenter));
   
    

    TimeMouseInCenter=timesMouseInCenter*FrameStep; % in seconds
    TimeMouseOutOfCenter=timesMouseOutOfCenter*FrameStep; % in seconds
    
    Thigmotaxis_Percent_90min=100*TimeMouseOutOfCenter/(TimeMouseOutOfCenter+TimeMouseInCenter);
    TimePercent_CenterZone_90min=100-Thigmotaxis_Percent_90min;
    
end