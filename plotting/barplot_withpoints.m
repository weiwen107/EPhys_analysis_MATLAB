
%Labels on x axis
X1 = 1;
X2 = 2;
% X3 = 3;
% X4 = 4;
% X5 = 5;
% X6 = 6;

%scale (for converting s into ms)
scale = 1000; 

%values for each label
Y1 = DR_CT_CNO_L.decaytau*scale;
aveY1 = nanmean(Y1);
Y2 = temp_1*scale;
aveY2 = nanmean(Y2);
% Y3 = adult_DR_CNO.amp*scale;
% aveY3 = nanmean(Y3);
% Y4 = shk3_WT_CNO_DW.Vm*scale;
% aveY4 = nanmean(Y4);
% Y5 = shk3_WT_DR_CNO_DW.Vm*scale;
% aveY5 = nanmean(Y5);
% Y6 = shk3_DR_CNO_DW.Vm*scale;
% aveY6 = nanmean(Y6);

figure() %for 2-3 groups

% figure('Position',[500 100 700 500]); %for 5-6 groups
bar(X1,aveY1,'EdgeColor','k','FaceColor','w','LineWidth',1,'BarWidth',0.4)

% hold on
% bar(X2,aveY2,'EdgeColor',[0.98 0.5 0.45],'FaceColor',[0.98 0.5 0.45],'LineWidth',1,'BarWidth',0.4)
% bar(X3,aveY3,'EdgeColor',[0.27 0.51 0.71],'FaceColor',[0.27 0.51 0.71],'LineWidth',1,'BarWidth',0.4)
% bar(X4,aveY4,'EdgeColor','k','FaceAlpha',0,'LineWidth',1,'BarWidth',0.4)
% bar(X4,aveY4,'EdgeColor',[0.98 0.5 0.45],'FaceColor',[0.98 0.5 0.45],'LineWidth',1,'BarWidth',0.4)
% bar(X5,aveY5,'EdgeColor',[0.27 0.51 0.71],'FaceColor',[0.27 0.51 0.71],'LineWidth',1,'BarWidth',0.4)
% bar(X4,aveY4,'EdgeColor',[0.96 0.64 0.38],'FaceColor',[0.96 0.64 0.38],'LineWidth',1,'BarWidth',0.4)
% bar(X5,aveY5,'EdgeColor',[0.53 0.81 0.98],'FaceColor',[0.53 0.81 0.98],'LineWidth',1,'BarWidth',0.4)
% bar(X6,aveY6,'EdgeColor',[0.27 0.51 0.71],'FaceColor',[0.27 0.51 0.71],'LineWidth',1,'BarWidth',0.4)
hold on
bar(X2,aveY2,'EdgeColor',[0.41 0.41 0.41],'FaceColor',[0.41 0.41 0.41],'LineWidth',1,'BarWidth',0.4)

ax = gca;
ax.FontSize = 12;
ax.LineWidth = 2;
ax.YLabel.String = 'pA';
%ax.YLim = [0 10e-4];
%ax.YTick = [0 50 100];
ax.XTick = [1 2];% 4 5 6];
ax.XTickLabel = [];
ax.XTickLabel = {'CT+','CT-'};%,'CNO','DR+CNO'};
ax.TickLabelInterpreter = 'none';
ax.XAxisLocation = 'origin';
%ax.XTickLabel.FontSize = 14;
%ax.XTickLabelRotation = 60;
%legend('hM4Di+CNO','Scaled','ctrl+CNO','FontSize',12,'Location','southeast')

%%%errorbars
%errhigh = amp_avg+amp_std;
%errlow = amp_avg-amp_std;

Y1_std = nanstd(Y1);
Y1_sem = Y1_std/sqrt(count_non_nan(Y1));
er1 = errorbar(X1,aveY1,Y1_sem);
er1.Color = 'k';
er1.LineWidth = 2;
er1.CapSize = 12;
X_1 = repmat(X1,1,numel(Y1));

Y2_std = nanstd(Y2);
Y2_sem = Y2_std/sqrt(count_non_nan(Y2));
er2 = errorbar(X2,aveY2,Y2_sem);
er2.Color = 'k';
er2.LineWidth = 2;
er2.CapSize = 12;
X_2 = repmat(X2,1,numel(Y2));

% Y3_std = nanstd(Y3);
% Y3_sem = Y3_std/sqrt(count_non_nan(Y3));
% er3 = errorbar(X3,aveY3,Y3_sem);
% er3.Color = 'k';
% er3.LineWidth = 2;
% er3.CapSize = 12;
% X_3 = repmat(X3,1,numel(Y3));

% Y4_std = nanstd(Y4);
% Y4_sem = Y4_std/sqrt(count_non_nan(Y4));
% er4 = errorbar(X4,aveY4,Y4_sem);
% er4.Color = 'k';
% er4.LineWidth = 2;
% er4.CapSize = 12;
% X_4 = repmat(X4,1,numel(Y4));
% 
% Y5_std = std(Y5);
% Y5_sem = Y5_std/sqrt(count_non_nan(Y5));
% er5 = errorbar(X5,aveY5,Y5_sem);
% er5.Color = 'k';
% er5.LineWidth = 2;
% er5.CapSize = 12;
% X_5 = repmat(X5,1,numel(Y5));
% 
% Y6_std = std(Y6);
% Y6_sem = Y6_std/sqrt(count_non_nan(Y6));
% er6 = errorbar(X6,aveY6,Y6_sem);
% er6.Color = 'k';
% er6.LineWidth = 2;
% er6.CapSize = 12;
% X_6 = repmat(X6,1,numel(Y6));

%generating x values for swarmplot
dp_offset = 0.2;
[X_1_s,Y1_s] = swarmplot(X_1+dp_offset,Y1,0.1);
[X_2_s,Y2_s] = swarmplot(X_2+dp_offset,Y2,0.1);
% [X_3_s,Y3_s] = swarmplot(X_3+dp_offset,Y3,0.1);
% [X_4_s,Y4_s] = swarmplot(X_4+dp_offset,Y4,0.1);
% [X_5_s,Y5_s] = swarmplot(X_5+dp_offset,Y5,0.1);
% [X_6_s,Y6_s] = swarmplot(X_6+dp_offset,Y6,0.1);

%scatter plot

scatter(X_1_s,Y1_s,'MarkerEdgeColor','k','MarkerFaceColor','k',...
    'LineWidth',1,'SizeData',20)
scatter(X_2_s,Y2_s,'MarkerEdgeColor',[0.8 0.8 0.8],'MarkerFaceColor',[0.8 0.8 0.8],...
    'LineWidth',1,'SizeData',20)
% scatter(X_3_s,Y3_s,'MarkerEdgeColor',[0.4 0.4 0.4],'MarkerFaceColor',[0.4 0.4 0.4],...
%     'LineWidth',1,'SizeData',20)
% 
% scatter(X_4_s,Y4_s,'MarkerEdgeColor','k','MarkerFaceColor','k',...
%     'LineWidth',1,'SizeData',20)
% scatter(X_4_s,Y4_s,'MarkerEdgeColor',[0.8 0.8 0.8],'MarkerFaceColor',[0.8 0.8 0.8],...
%     'LineWidth',1,'SizeData',20)
% scatter(X_5_s,Y5_s,'MarkerEdgeColor',[0.4 0.4 0.4],'MarkerFaceColor',[0.4 0.4 0.4],...
%     'LineWidth',1,'SizeData',20)
% scatter(X_6_s,Y6_s,'MarkerEdgeColor',[0.4 0.4 0.4],'MarkerFaceColor',[0.4 0.4 0.4],...
%     'LineWidth',1,'SizeData',25)
% scatter(X_4_s,Y4_s,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],...
%     'LineWidth',1,'SizeData',25)
% scatter(X_5_s,Y5_s,'MarkerEdgeColor',[0.2 0.2 0.2],'MarkerFaceColor',[0.2 0.2 0.2],...
%     'LineWidth',1,'SizeData',25)
title('New CT data')
hold off
box off