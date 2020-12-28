%%%%
%This script imports the waveform average data of each cell and take their
%average to get metaaverage waveform (for each condition)

%% import the waveform average of each cell, for each condition, and calculate the meta average
save_results = 0;

figure_on_unscaled = 1;

figure_on_scaled = 1;

fp_analyzed_data = 'C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di';

import_file_name = 'waveform_average_rise1.xlsx';

save_filename = 'chronic_DREADDs_waveform_average_shk3_wt_dw_rise1.mat';

% experimental conditions (should corresponds to sheet names in the excel
% file)
condn = {'shk3_wt_CNO_DW', 'shk3_wt_DR_CNO_DW', 'shk3_CNO_DW', 'shk3_DR_CNO_DW'};

fn = condn;

data_temp = NaN(150,30);

waveform_ave = struct;

for fi = 1 : numel(fn)
    waveform_ave.(fn{fi}) =...,
        xlsread(strcat(fp_analyzed_data,'\',import_file_name), fn{fi});
end

%calculate average

for fii = 1 : numel(fn)
    for di = 1:size(waveform_ave.(fn{fii}),1)
        meta_ave.(fn{fii})(di,1) = nanmean(waveform_ave.(fn{fii})(di,:));
    end
end

%% Scale peaks: scale condition2 to condition1

condition1 = meta_ave.shk3_CNO_DW;
condition2 = meta_ave.shk3_DR_CNO_DW;

[peak_cond1, peak_ind_cond1] = min(condition1);
[peak_cond2, peak_ind_cond2] = min(condition2);

scale_factor = mean(condition1((peak_ind_cond1-1):(peak_ind_cond1+1),1))/...,
    mean(condition2((peak_ind_cond2-1):(peak_ind_cond2+1),1));

condition2_scaled = NaN(size(condition1,1),1);
diff = size(condition1,1)-size(condition2,1);

if diff >=0
    condition2_scaled(diff+1:end) = condition2 .* scale_factor;
else
    condition2_scaled = condition2(abs(diff)+1 : end) .* scale_factor;
end
    

%plotting unscaled
if figure_on_unscaled == 1
    figure();
    plot(condition1(10:110),'Color',[0.98 0.5 0.45],'LineWidth',4)
    hold on
    plot(condition2(10:110),'Color',[0.27 0.51 0.71],'LineWidth',4)
    %hold on
    %plot(meta_ave.DR_CNO_24(20:110),'Color',[0.27 0.51 0.71],'LineWidth',4)
    title('unscaled')
    box off
    hold off
end

    %plotting scaled
if figure_on_scaled == 1    
    figure();
    plot(condition1(10:110),'Color',[0.98 0.5 0.45],'LineWidth',4)
    hold on
    plot(condition2_scaled(10:110),'Color',[0.27 0.51 0.71],'LineWidth',4)
%     hold on
%     plot(shk3_scaled(20:110),'Color',[0.27 0.51 0.71],'LineWidth',4)
    title('scaled')
    %legend('WT','shk3');
    box off
    hold off
end

%% save data
if save_results == 1
    cd (fp_analyzed_data)
    save(save_filename,'waveform_ave','meta_ave')
end