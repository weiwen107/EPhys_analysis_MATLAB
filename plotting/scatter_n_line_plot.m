table_num = 9;

%choose which property to plot
curr_ta = decaytau;
curr_ta.Variables = curr_ta.Variables .*1000;

%set x labels
X = 1:2;%(1:size(curr_ta,2));

%individual cell status
cst = NaN(size(curr_ta,1),1);

%% plotting
figure();
for xi = 1:size(curr_ta,1)
    %Rin_change = abs(diff(Rin(xi,1:2).Variables))/nanmean(Rin(xi,1:2).Variables);
    Ra_change = abs(diff(Ra(xi,1:2).Variables))/nanmean(Ra(xi,1:2).Variables);
    amp_diff = diff(amp(xi,1:2).Variables);
    if  Ra_change <= 0.2 %&& amp_diff <=0
        cst(xi) = xi;
        plot(X,curr_ta(xi,1:2).Variables,'-o','Color','k','MarkerSize',10,...
            'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',1)
        hold on
    end
end
%hold off
%box off
%format the plot
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 2;
ax.YLabel.String = 'decay (ms)';
ax.YLabel.FontSize = 14;
ax.YMinorTick = 'off';
%ax.YLim = [0 15];
%ax.YTick = [0 5 10 15];
% ax.YTickLabel = [0 5 10 15];
ax.XLim = [0 3];
ax.XTick = [1 2];% 3];
ax.XTickLabel = {'ACSF', 'Naspm'};% 'After'};
ax.TickLabelInterpreter = 'none';
%ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

hold off
box off

%% test
sample1 = NaN(size(curr_ta,1),1);
sample2 = NaN(size(curr_ta,1),1);


for tii = 1:size(cst,1)
    if isnan(cst(tii))
        sample1(tii) = NaN;
        sample2(tii) = NaN;
    else
        sample1(tii) = curr_ta(tii,1).Variables;
        sample2(tii) = curr_ta(tii,2).Variables;
    end
end

[p,h,stats] = ranksum(sample1, sample2);

