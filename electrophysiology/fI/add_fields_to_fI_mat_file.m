%fI_grouped_data_name
saved_file_name = 'chronic_hm4di_24_5_fI.mat';

%location of grouped files
fp_grouped_data = ...,
    'C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di\fI_data_by_groups';

%experimental conditions
exp_con = {'WT','WT_CNO','DR_CNO','WT_saline'};

%confirm the condition of each column first!
fp_import_file = 'C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di\matlab_import_file';
import_filename = 'matlab_import_fI.xlsx';

cd(fp_grouped_data)
load(saved_file_name)

%%
fields = {'Vm'};
data_temp = cell(1,numel(fields));
struct_temp = cell(1,numel(exp_con));



for shi = 1:numel(fields)
    data_temp{shi} = xlsread(strcat(fp_import_file,'\',import_filename),fields{shi});
end

for fi = 1:numel(fields)
    
    for sti = 1:size(struct_temp,2)
                 
        if sti == 1
            WT.(fields{fi}) = data_temp{fi}(1:end,sti);
        elseif sti == 2
            WT_CNO.(fields{fi}) = data_temp{fi}(1:end,sti);
        elseif sti == 3
            DR_CNO.(fields{fi}) = data_temp{fi}(1:end,sti);
        elseif sti == 4
            WT_saline.(fields{fi}) = data_temp{fi}(1:end,sti);
        end

    end   
end

%% 

save(saved_file_name,exp_con{1},exp_con{2},exp_con{3},exp_con{4},'current_inj')