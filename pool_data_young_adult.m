% This script pools data from separate experimental conditions 
% into a single group for the purpose of comparison  
% individual experimental conditions will be saved in 'sample'
% pooled data will be saved in 'pooled_sample'

%% 

save_results = 1;

save_file_name = 'mini_CP_adult_NT_CNO_pooled.mat';

save_file_location = ...,
    'C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di\mini_data_by_groups';

%% data extraction and grouping
%name of groups
%the first level indicates name of the experimental group (juvenile/adult)
%the second level indicates the experimental conditions under each group
%that will be combined into a single column (NT/CNO)
groups = {{'WT','WT_CNO'}, {'adult_WT','adult_CNO'}};

%name of the variable 
parameter = {'amp','frq','risetime','decaytau','charge','Vm','Ra','Cm','Rin'};

%testing mode: 0 for single column variable, 1 for multiple column
%variable
test_mode = 0;

%data selection mode: only applicable when plot_mode == 1
%1 for no normalization(start from the 1st current step)
%2 for normalization (start from the rheobase step)
select_mode = 2;

%stimulation amplitude: only applicable when plot_mode == 1
%indicates the stimulation amplitude (current injected, in pA)
stim = 0;

%confidence level
cl = 0.05;

sample = cell(1,numel(groups));
pooled_sample = cell(1,numel(groups));

for pi = 1:numel(parameter) 
    for gi = 1:numel(groups)
        for eci = 1:numel(groups{1,gi})
            if test_mode == 0
                    sample{1,gi}{1,eci} = eval(strcat(groups{1,gi}{1,eci},'.',parameter{pi}));
                elseif test_mode ==1
                    data_temp = eval(strcat(groups{1,gi}{1,eci},'.',parameter{pi}));
                    rheo = eval(strcat(groups{1,gi{1,eci}},'.','Rheobase'));
                    sample{1,gi}{1,eci} = get_step_values(data_temp, rheo, select_mode, stim);
            end

            if eci == 1
                pooled_sample{1,gi}.(parameter{pi}) = sample{1,gi}{1,1};
            else
                pooled_sample{1,gi}.(parameter{pi}) = ...,
                    cat(1, pooled_sample{1,gi}.(parameter{pi}), sample{1,gi}{1,eci});
            end
        end
    end
end

CP_pooled = pooled_sample{1,1};
adult_pooled = pooled_sample{1,2};

%% save
if save_results == 1
    cd(save_file_location)
    save(save_file_name, 'CP_pooled','adult_pooled')
end

