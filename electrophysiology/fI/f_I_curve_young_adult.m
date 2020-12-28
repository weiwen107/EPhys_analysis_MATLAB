%% f-I curve

%plot on
plot_on = 1;

% data range for plotting
plot_range = 1:20;
plot_variable = 'mean_IFR';
plot_stat = 'ave_y_a';
err = 'sem_y_a';

curr_y_y = eval(strcat(plot_stat,'.young.',plot_variable));
curr_y_err_y = eval(strcat(err,'.young.',plot_variable));
curr_y_a = eval(strcat(plot_stat,'.adult.',plot_variable));
curr_y_err_a = eval(strcat(err,'.adult.',plot_variable));

if plot_on == 1

    figure
    hold on
    
    %young
    errorbar(current_inj(plot_range),curr_y_y(plot_range,3),curr_y_err_y(plot_range,3),'-o','MarkerSize',8,...
            'MarkerEdgeColor','k', 'MarkerFaceColor','k','LineWidth',2,...
            'Color','k')
    %adult
    errorbar(current_inj(plot_range),curr_y_a(plot_range,1),curr_y_err_a(plot_range,1),'-o','MarkerSize',8,...
            'MarkerEdgeColor',[105 105 105]./255, 'MarkerFaceColor',[105 105 105]./255,'LineWidth',2,...
            'Color',[105 105 105]./255)


    ax = gca;
    ax.FontSize = 12;
    ax.LineWidth = 2;
    ax.YLabel.String = 'mean IFR (Hz)';
    ax.YLabel.FontSize = 14;
    ax.XLabel.String = 'Injected Current (pA)';
    ax.YLabel.FontSize = 14;
    ax.YLabel.Interpreter = 'none';
    ax.XLim = [0 450];
    ax.YLim = [0 inf];

    legend('P24-29','P45-50','FontSize',12,'Location','northwest','Interpreter','none','Box','off')

end

%% first AP trace overlay
condition1 = meta_ap_average.young{1,2};
condition2 = meta_ap_average.adult{1,2};

plot_range = (1:80);

if plot_on == 1
    figure();
    plot(condition1(plot_range),'Color','k','LineWidth',2)
    hold on
    plot(condition2(plot_range),'Color',[105 105 105]./255,'LineWidth',2)
    %hold on
    %plot(meta_ave.DR_CNO_24(20:110),'Color',[0.27 0.51 0.71],'LineWidth',4)
    title('First AP at rheobase step')
    legend('P24-29','P45-50','FontSize',12,'Location','northeast','Interpreter','none','Box','off')
    box off
    hold off
end