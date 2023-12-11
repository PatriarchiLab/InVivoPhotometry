%Bar plots with the individual data points for different cohorts
%You can use Prism Program instead :) 

%author: Chernysheva Maria
%version:01.09.22

function SEM_Input = bar_plots_cohorts(input_data,legend_cohort,color_cohort, plot_title,YLabel,font_size, paired)

%some arbitrary numbers for the plot's design:
 
%Figure Proportions
pbaspect_1 = 1;
pbaspect_2 = 1.5;
pbaspect_3 = 1;
step2axis = -0.3;

data_marker_size = font_size*0.75; %how big are the individual grey dots
error_bar_width = 6;
error_bar_line_idth = 1.5;
bar_width = 0.35;
line_width = 0.5;
color_grey = [0.9 0.9 0.9];

for cohort_number = 1:size(input_data,2)    
SEM_Input(cohort_number) = std(input_data{cohort_number})/sqrt(length(input_data{cohort_number}));
x(cohort_number) = step2axis*cohort_number+0.8*cohort_number;
end


for cohort_number = 1:size(input_data,2)
bar(x(cohort_number),mean(input_data{cohort_number}),'FaceColor',color_cohort{cohort_number},'EdgeColor','none','BarWidth', bar_width);
hold on
e=errorbar(x(cohort_number),mean(input_data{cohort_number}),SEM_Input(cohort_number),'LineWidth',error_bar_line_idth);
e.Color = color_cohort{cohort_number};
e.CapSize=error_bar_width;
hold on
if paired==0
    
for cohort_number = 2:size(input_data,2)
% Calculate P-value for columns, it gives the quick overview, but better use
% PRISM to quantify it properly (with multiple comparisons)
[p_value(cohort_number-1),h]=ranksum (input_data{1},input_data{cohort_number});
string_p_value{cohort_number-1}=['p_1_-_',num2str(cohort_number),'=',num2str(round(p_value(cohort_number-1),4))];
end
for cohort_number = 1:size(input_data,2)
%plot circles for the individual data points
scatter (ones(1,length(input_data{cohort_number}))*x(cohort_number),input_data{cohort_number},data_marker_size,'MarkerEdgeColor', color_cohort{cohort_number},'MarkerFaceColor','w')
hold on
end
end
end  

if paired == 1
    for cohort_number=1:size(input_data,2)-1
    
% Calculate P-value for columns, it gives the quick overview, but better use
% PRISM to quantify it properly (with multiple comparisons)
[p_value(cohort_number),h]=ranksum (input_data{1},input_data{cohort_number+1});
% string_p_value{cohort_number}=['p_1_-_',num2str(cohort_number+1),'=',num2str(round(p_value(cohort_number),4)),' '];
string_p_value{cohort_number}=['p',' = ',num2str(round(p_value,4)),' '];

    for Current_mouse_InputData=1:length(input_data{1})
    %plot dots representing "mice"
    plot ([x(cohort_number) x(cohort_number+1)],[input_data{cohort_number}(Current_mouse_InputData) input_data{cohort_number+1}(Current_mouse_InputData)],'.k','LineWidth',line_width,'Color',color_grey,'MarkerSize',data_marker_size,'MarkerFaceColor',color_grey)
    hold on
    %plot the line connecting the "mice" dots
    line ([x(cohort_number) x(cohort_number+1)],[input_data{cohort_number}(Current_mouse_InputData) input_data{cohort_number+1}(Current_mouse_InputData)],'LineWidth',line_width,'Color',color_grey)
    hold on
    end
    end
end

xlim([0 (x(end)+x(1))])
set(gca,'XTick',[x(:)]);


title ([plot_title,{string_p_value{:}},{''}],'FontName','Arial','FontWeight','light','Fontsize',font_size,'Interpreter','none')
% title ([plot_title,{''}],'FontName','Arial','FontWeight','light','Fontsize',font_size-2)
ylabel(YLabel,'Fontsize',font_size,'Color','k','FontName','Arial');
set(gca,'TickDir','out');
set(gca,'TickLength',[0.02, 0.08])
xticklabels({legend_cohort{:}})
xtickangle(45)

pbaspect([pbaspect_1 pbaspect_2 pbaspect_3])
box off

end