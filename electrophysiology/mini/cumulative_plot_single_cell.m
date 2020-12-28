%% Choose cells from both conditions
%whether to save results
save_file = 1;

%name of the cumulative file
save_file_name = 'adult_cumulative_data';

%location of mini analysis files
fp_analyzed_data = ...,
    'C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di\analyzed_mini_results\rise_1\';

%experiments where each cells are from
mini_exp = load(strcat(fp_analyzed_data, 'MINIANALYSIS_191115.mat'));
AMP_ALL_exp = mini_exp.AMP_ALL;

mini_ctrl = load(strcat(fp_analyzed_data, 'MINIANALYSIS_191115.mat'));
AMP_ALL_ctrl = mini_ctrl.AMP_ALL;

% cell number within the corresponding experiment
cell_num = [1 5];

%% data preparation
exp_con = NaN(size(AMP_ALL_exp{1,1}{1,cell_num(1)}));
ctrl_con = NaN(size(AMP_ALL_ctrl{1,1}{1,cell_num(2)}));

tracect_exp = 1;
tracect_ctrl = 1;

for eci = 1:size(AMP_ALL_exp{1,1}{1,cell_num(1)},2)
    if sum(~isnan(AMP_ALL_exp{1,1}{1,cell_num(1)}(:,eci))) == 0 || ...
            nanmean(AMP_ALL_exp{1,1}{1,cell_num(1)}(:,eci)) == 0
        continue
    else       
        lgth = size(AMP_ALL_exp{1,1}{1,cell_num(1)}(:,eci),1);
        exp_con(1:lgth,tracect_exp) = AMP_ALL_exp{1,1}{1,cell_num(1)}(:,eci);
        tracect_exp = tracect_exp + 1;
    end
end

for cci = 1:size(AMP_ALL_ctrl{1,1}{1,cell_num(2)},2)
    if sum(~isnan(AMP_ALL_ctrl{1,1}{1,cell_num(2)}(:,cci))) == 0 || ...
            nanmean(AMP_ALL_ctrl{1,1}{1,cell_num(2)}(:,cci)) == 0
        continue
    else
        lgth = size(AMP_ALL_ctrl{1,1}{1,cell_num(2)}(:,cci),1);
        ctrl_con(1:lgth,tracect_ctrl) = AMP_ALL_ctrl{1,1}{1,cell_num(2)}(:,cci);
        tracect_ctrl = tracect_ctrl + 1;
    end
end

%convert NaNs into zeros, then remove zeros and concatenate into a single column.

for ei = 1:size(exp_con,2)
    for ej = 1:size(exp_con,1)
        if isnan(exp_con(ej,ei))
            exp_con(ej,ei) = 0;
        end
    end
end

exp_con_nonzero = nonzeros(exp_con);

for ci = 1:size(ctrl_con,2)
    for cj = 1:size(ctrl_con,1)
        if isnan(ctrl_con(cj,ci))
            ctrl_con(cj,ci) = 0;
        end
    end
end

ctrl_con_nonzero = nonzeros(ctrl_con);

%% pick events and calculate cumulative distribution
%randomly pick 300 events from each condition
exp_con_final = exp_con_nonzero(1:300);

% rand_index_exp = randi([1 size(exp_con_nonzero,1)],200,1);
% exp_con_final = exp_con_nonzero(rand_index_exp);

ctrl_con_final = ctrl_con_nonzero(1:300);

% rand_index_ctrl = randi([1 size(ctrl_con_nonzero,1)],300,1);
% ctrl_con_final = ctrl_con_nonzero(rand_index_ctrl);

%generate cumulative statistics
min_exp = min(exp_con_final);
max_exp = max(exp_con_final);
[exp_X, exp_Y] = cumhist(exp_con_final,[min_exp max_exp],0.01);
min_ctrl = min(ctrl_con_final);
max_ctrl = max(ctrl_con_final);
[ctrl_X, ctrl_Y] = cumhist(ctrl_con_final,[min_ctrl max_ctrl],0.01);

exp_sort = sort(exp_con_final);
ctrl_sort = sort(ctrl_con_final);

%use linear fitting to get scale factor
f = fittype('a*x+b','coefficients',{'a','b'});
[lnfit,gof] = fit(ctrl_sort,exp_sort,f,'StartPoint',[min_exp min_ctrl]);
exp_scaled = (exp_sort-lnfit.b) ./ lnfit.a;

%plot scaled data
figure()
plot(ctrl_sort,ctrl_sort,'k', 'Linewidth',2)
hold on
plot(ctrl_sort,exp_sort,'o','MarkerSize',8)
plot(lnfit)

ax1 = gca;
ax1.FontSize = 12;
ax1.LineWidth = 2;
ax1.YLabel.String = 'DR+CNO (pA)';
ax1.YLabel.FontSize = 14;
ax1.XTick = [0 20 40 60 80];
ax1.XLabel.String = 'CNO Only (pA)';
ax1.XLabel.FontSize = 14;
legend('CNO vs. CNO','DR+CNO vs. CNO','best fit','FontSize',12,'Location','southeast')
legend('boxoff')
title('P45-50')
box off

hold off

%average the scale factor for each pair of corresponding points (sorted
%vector)
% scale_factor = mean(exp_con_final./ctrl_con_final);
% exp_scaled = exp_sort ./ scale_factor;

%two sample Kolmogorov-Smirnov test
[exp_X_scaled, exp_Y_scaled] = cumhist(exp_scaled,[min(exp_scaled) max(exp_scaled)],0.01);
[h1,p1,ks2stat1] = kstest2(exp_scaled,ctrl_sort);
[h2,p2,ks2stat2] = kstest2(exp_sort,ctrl_sort);

%% plotting cumulative plot
figure();
plot(exp_X, smooth(exp_Y),'Color',[0.27 0.51 0.71],'LineWidth',2)
hold on
% exp_X_scaled = cat(2,exp_X_scaled, max_exp);
% exp_Y_scaled = cat(2,exp_Y_scaled,100);
plot(exp_X_scaled, smooth(exp_Y_scaled),'k','LineWidth',2)

% ctrl_X = cat(2,ctrl_X,max_exp);
% ctrl_Y = cat(2,ctrl_Y,100);
plot(ctrl_X, smooth(ctrl_Y),'Color',[0.98 0.5 0.45],'LineWidth',2)

ax = gca;
ax.FontSize = 12;
ax.LineWidth = 2;
ax.YLabel.String = 'Cumulative %';
ax.YLabel.FontSize = 14;
ax.XLim = [0 50];
ax.YLim = [0 110];
ax.YTick = [0 50 100];
ax.XLabel.String = 'mEPSC amplitude (pA)';
ax.XLabel.FontSize = 14;
legend('DR+CNO','Scaled','CNO Only','FontSize',12,'Location','southeast')
legend('boxoff')
box off
hold off

%% save to file
if save_file == 1
    cd('C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di\cumulative_data')
    save(save_file_name,'exp_X','exp_Y','ctrl_X','ctrl_Y','exp_X_scaled','exp_Y_scaled',...
        'h1','p1','ks2stat1','h2','p2','ks2stat2')
end