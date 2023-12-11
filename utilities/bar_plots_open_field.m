%Bar plots with the individual data points for different cohorts
%You can use Prism Program instead :) 

%author: Chernysheva Maria
%version:01.09.22

function SEM_Input = bar_plots_open_field(InputData,Legend_Cohort,Color_Cohort, Title,YLabel,FontSize, paired)

%some arbitrary numbers for the plot's design:
 
%Figure Proportions
pbaspect_1=1.2;
pbaspect_2=2;
pbaspect_3=1;
step2axis=-0.3;

DataMarkerSize=10; %how big are the individual grey dots
ErrorBarWidth=6;
ErrorBarLineWidth=1.5;
barWidth=0.35;
LineWidth=0.5;
Color_Grey=[0.9 0.9 0.9];

for Cohort_Number=1:size(InputData,2)    
SEM_Input(Cohort_Number)=std(InputData{Cohort_Number})/sqrt(length(InputData{Cohort_Number}));
x(Cohort_Number)=step2axis*Cohort_Number+0.8*Cohort_Number;


end


for Cohort_Number=1:size(InputData,2)
bar(x(Cohort_Number),mean(InputData{Cohort_Number}),'FaceColor',Color_Cohort{Cohort_Number},'EdgeColor','none','BarWidth', barWidth);
hold on
e=errorbar(x(Cohort_Number),mean(InputData{Cohort_Number}),SEM_Input(Cohort_Number),'LineWidth',ErrorBarLineWidth);
e.Color = Color_Cohort{Cohort_Number};
e.CapSize=ErrorBarWidth;
hold on
if paired==0
    
for Cohort_Number=2:size(InputData,2)
% Calculate P-value for columns, it gives the quick overview, but better use
% PRISM to quantify it properly (with multiple comparisons)
[p_value(Cohort_Number-1),h]=ranksum (InputData{1},InputData{Cohort_Number});
string_p_value{Cohort_Number-1}=['p_1_-_',num2str(Cohort_Number),'=',num2str(round(p_value(Cohort_Number-1),4))];
end
for Cohort_Number=1:size(InputData,2)
%plot circles for the individual data points
scatter (ones(1,length(InputData{Cohort_Number}))*x(Cohort_Number),InputData{Cohort_Number},DataMarkerSize,'MarkerEdgeColor', Color_Cohort{Cohort_Number},'MarkerFaceColor','w')
hold on
end
end
end  

if paired==1
    for Cohort_Number=1:size(InputData,2)-1
    
% Calculate P-value for columns, it gives the quick overview, but better use
% PRISM to quantify it properly (with multiple comparisons)
% [p_value(Cohort_Number),h]=signrank (InputData{1},InputData{Cohort_Number+1});
% string_p_value{Cohort_Number}=['p_1_-_',num2str(Cohort_Number+1),'=',num2str(round(p_value(Cohort_Number),4)),' '];

    for Current_mouse_InputData=1:length(InputData{1})
    %plot dots representing "mice"
    plot ([x(Cohort_Number) x(Cohort_Number+1)],[InputData{Cohort_Number}(Current_mouse_InputData) InputData{Cohort_Number+1}(Current_mouse_InputData)],'.k','LineWidth',LineWidth,'Color',Color_Grey,'MarkerSize',DataMarkerSize,'MarkerFaceColor',Color_Grey)
    hold on
    %plot the line connecting the "mice" dots
    line ([x(Cohort_Number) x(Cohort_Number+1)],[InputData{Cohort_Number}(Current_mouse_InputData) InputData{Cohort_Number+1}(Current_mouse_InputData)],'LineWidth',LineWidth,'Color',Color_Grey)
    hold on
    end
    end
end

xlim([0 (x(end)+x(1))])
set(gca,'XTick',[x(:)]);


% title ([Title,{string_p_value{:}},{''}],'FontName','Arial','FontWeight','light','Fontsize',FontSize-2)
title ([Title,{''}],'FontName','Arial','FontWeight','light','Fontsize',FontSize-2)
ylabel(YLabel,'Fontsize',FontSize,'Color','k','FontName','Arial');
set(gca,'TickDir','out');
set(gca,'TickLength',[0.025, 0.1])
xticklabels({Legend_Cohort{:}})
xtickangle(45)

pbaspect([pbaspect_1 pbaspect_2 pbaspect_3])
box off

end