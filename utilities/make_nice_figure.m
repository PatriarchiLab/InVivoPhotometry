function make_nice_figure(name_plot,fontSize_4figures)

box off
set(gca,'TickDir','out');
set(gca,'FontSize',fontSize_4figures)
set(gca, 'Color', 'None');
save_plots(name_plot)

end