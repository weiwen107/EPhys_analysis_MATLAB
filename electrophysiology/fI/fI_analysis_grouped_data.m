%%%% This script takes grouped fI data sets and conducts group analysis

%% Group info
%location where the mat file will be saved
fp_analyzed_data = 'C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di\fI_data_by_groups';

%name of the mat file
filename = 'hM4Di_shk3_grouped_analysis.mat';

%save results
save_results = 1;

%name of groups (experimental conditions)
cond = {'shk3_CNO','shk3_DR_CNO'};

%% import grouped fI data 
%groups are sequentially stored as indicated in cond
group_struct = cell(1,numel(cond));

for gi = 1:numel(cond)
    group_struct{1,gi} = eval(cond{1,gi});
end

%% Bursting neurons
% Bursts are defined as two spikes having an interspike interval less than
% 10 ms.

% quantifing IFR at 300pA current injection because: 
% 1. at 200pA, almost all cells have IFRs lower than 100 Hz
% 2. at 400 pA, almost all cells have IFRs higher than 100 Hz
% 300 pA seems like a middle ground
% same reason for choosing 400pA for mean IFR

stim_IFR = 300;
stim_mIFR = 400;

IFR_burst_cell_counter = zeros(1,numel(cond));
IFR_burst_cell_index = cell(1,numel(cond)); %binary indicator: 0 as non-burst, 1 as burst
IFR_prct_burst = zeros(1,numel(cond));

mIFR_burst_cell_counter = zeros(1,numel(cond));
mIFR_burst_cell_index = cell(1,numel(cond)); %binary indicator: 0 as non-burst, 1 as burst
mIFR_prct_burst = zeros(1,numel(cond));

%IFR
for gii = 1:numel(cond)
    for ci = 1:size(group_struct{1,gii}.IFR,2)
        if isnan(group_struct{1,gii}.IFR(1,ci))
            continue
        else
            curr_IFR = group_struct{1,gii}.IFR(stim_IFR/20,ci);
            if curr_IFR > 100 || curr_IFR == 100
                IFR_burst_cell_counter(1,gii) = IFR_burst_cell_counter(1,gii) + 1;
                IFR_burst_cell_index{1,gii}(1,ci) = 1;
            else
                IFR_burst_cell_index{1,gii}(1,ci) = 0;
            end
        end
    end
    
    IFR_prct_burst(1,gii) = IFR_burst_cell_counter(1,gii) / size(IFR_burst_cell_index{1,gii},2) * 100;
end

%mean IFR
for gii = 1:numel(cond)
    for ci = 1:size(group_struct{1,gii}.mean_IFR,2)
        if isnan(group_struct{1,gii}.mean_IFR(1,ci))
            continue
        else
            curr_mIFR = group_struct{1,gii}.mean_IFR(stim_mIFR/20,ci);
            if curr_mIFR > 100 || curr_mIFR == 100
                mIFR_burst_cell_counter(1,gii) = mIFR_burst_cell_counter(1,gii) + 1;
                mIFR_burst_cell_index{1,gii}(1,ci) = 1;
            else
                mIFR_burst_cell_index{1,gii}(1,ci) = 0;
            end
        end
    end
    
    mIFR_prct_burst(1,gii) = mIFR_burst_cell_counter(1,gii) / size(mIFR_burst_cell_index{1,gii},2) * 100;
end

%% Correlation between fI curve slope and input resistance (using mean IFR)

fI_slope = cell(1,numel(cond));

%r- correlation coefficients
%p- corresponding p values
corr_val = cell(1,numel(cond));

figure_on = 1;

for gii = 1:numel(cond)
    for ci = 1:size(group_struct{1,gii}.mean_IFR,2)
        if isnan(group_struct{1,gii}.mean_IFR(1,ci))
            continue
        else
            curr_mIFR = group_struct{1,gii}.mean_IFR(1:20,ci);
            fit_start = find(curr_mIFR,1,'first');
            curr_X = (fit_start*20:20:400)';
            curr_Y = curr_mIFR(fit_start:20);
            
            curr_fit = fitlm(curr_X,curr_Y);
            fI_slope{1,gii}(ci,1) = curr_fit.Coefficients{2,'Estimate'};
        end
    end
    
    curr_Rin = group_struct{1,gii}.Rin(1:size(fI_slope{1,gii},1));
    [r,p] = corrcoef(curr_Rin,fI_slope{1,gii});
    corr_val{1,gii}.r = r;
    corr_val{1,gii}.p = p;
    
    if figure_on == 1
        figure()
        scatter(curr_Rin,fI_slope{1,gii})
        hold on
        fitt = fitlm(curr_Rin,fI_slope{1,gii});
        plot(fitt)
        hold off
    end
end


%% save results

if save_results == 1

    cd (fp_analyzed_data)
    save(filename,'stim_IFR','stim_mIFR','IFR_burst_cell_counter','IFR_burst_cell_index','IFR_prct_burst',...
        'mIFR_burst_cell_counter','mIFR_burst_cell_index','mIFR_prct_burst',...
        'fI_slope','corr_val')
end

