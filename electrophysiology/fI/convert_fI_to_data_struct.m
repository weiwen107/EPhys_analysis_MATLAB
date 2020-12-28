% This script loops through the analyzed_fI_results folder and groups
% them by experimental conditions: 
% For example: (these indice should match 'Cat' in the cell_id excel file) 
    % DR+CNO- 1
    % CNO- 2
    % NT- 3
    % NT+saline- 4

% The experimental condition of a cell can be uniquely determined by the date, 
% and the cell_num assigned on that day via a look-up cell_id_index table.  
%% %% Change these accordingly based on how you want to group the data
% in each anaylzed_fI_results file, dates can be extracted from the last six digits

%save results
save_results = 1;

%experiment name (correpsonding to the sheet name in the cell_id_index
% excel file)
exp_name = 'WT_adult_DW_48CNO';

%file name
saved_file_name = 'chronic_hm4di_fI_DW_adult_48h.mat';

%where to save grouped files
fp_grouped_data = ...,
    'C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di\fI_data_by_groups';

%experimental conditions
exp_con = {'adult_DR_CNO_48h','adult_CNO_48h'};

%import cell_id_index table 
cd('C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di\fI_data_by_groups')
cell_id_index = readtable('cell_id_index.xlsx','Sheet',exp_name);

%location of analyzed fI results
fp_analyzed_fI = ...,
    'C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di\analyzed_fI_results';

%% extract and group data
data_struct_temp = cell(1,numel(exp_con));

cd(fp_analyzed_fI)
all_files = dir;
file_num = numel(all_files)-2;

for fi = 1:file_num
    curr_name = all_files(fi+2).name;
    if ismember('m',curr_name(4:9))
        continue
    else
        curr_date = curr_name(4:9);
    end
    
    
    cy = strcat('20',curr_date(1:2));
    cM = curr_date(3:4);
    cday = curr_date(5:6);
    
    date = datetime(strcat(cy,'-',cM,'-',cday),'InputFormat','y-M-d',...
        'Format', 'M/d/y');
    if ~(ismember(date,cell_id_index.Date))
        continue
    else
        load(all_files(fi+2).name)
        
        for ci = 1:size(cell_id,2)
            if isempty(cell_id{1,ci})
                continue
            else
                if isempty(find(cell_id_index.Date == date & cell_id_index.Cell_num == ci,1))
                    continue
                else
                
                    row = find(cell_id_index.Date == date & cell_id_index.Cell_num == ci);
                    cond_i = cell_id_index{row,'Cat'};
                    cell_ID = cell_id_index{row,'Cell_ID'};
                    ti_start = cell_id{1,ci}(1,1);
                    ti_end = cell_id{1,ci}(20,1);
                    
                    data_struct_temp{cond_i}.MFR(1:20,cell_ID) = MFR{1,ci}(ti_start:ti_end,1);                 
                    data_struct_temp{cond_i}.IFR(1:20,cell_ID) = IFR{1,ci}(ti_start:ti_end,1);
                    data_struct_temp{cond_i}.mean_IFR(1:20,cell_ID) = IFR_ave{1,ci}(ti_start:ti_end,1);
                    data_struct_temp{cond_i}.Threshold(1:20,cell_ID) = V_th_1st{1,ci}(ti_start:ti_end,1);
                    data_struct_temp{cond_i}.ADP_ind(1:20,cell_ID) = adp_index{1,ci}(ti_start:ti_end,1);
                    data_struct_temp{cond_i}.Lat(1:20,cell_ID) = lat{1,ci}(ti_start:ti_end,1);
                    data_struct_temp{cond_i}.Udratio(1:20,cell_ID) = udratio_median{1,ci}(ti_start:ti_end,1);
                    data_struct_temp{cond_i}.Width_first(1:20,cell_ID) = width_1st{1,ci}(ti_start:ti_end,1);
                    data_struct_temp{cond_i}.Width_med(1:20,cell_ID) = width_median{1,ci}(ti_start:ti_end,1);
                    data_struct_temp{cond_i}.Rheobase(cell_ID,1) = rheobase(ci,1);
                end
            end
        end
    end
end
    
adult_DR_CNO_48h = data_struct_temp{1,1};
adult_CNO_48h = data_struct_temp{1,2};
current_inj = (20:20:400)';
%% save grouped data
if save_results == 1
    cd(fp_grouped_data)
    save(saved_file_name,exp_con{1},exp_con{2},'current_inj')
end
