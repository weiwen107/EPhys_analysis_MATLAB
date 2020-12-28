%% import from excel
%location where the import file will be saved
fp_import_file = 'C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di\matlab_import_file';
import_filename = 'matlab_import_fI_shk3.xlsx';

%location where the mat file will be saved
fp_analyzed_data = 'C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di\fI_data_by_groups';

%name of the mat file
filename = 'chronic_hm4di_fI_DW_shk3.mat';

%save results
save_results = 1;

%data location corresponding to each expermental condition
%the starting and end row # in excel file
%in order: DR_CNO, WT_CNO
dl = [2,21; 23,42];%; 44,63; 65,84];

%experimental conditions
% adult_WT_CNO = struct;
% adult_DR_CNO = struct;
% WT = struct;
% WT_saline = struct;

%field names (these should correspond to the sheet order in the import
%excel file, but the sheet names can be different)
fields = {'MFR', 'IFR', 'mean_IFR', 'Threshold', 'ADP_ind','Lat',...
    'Udratio','Width_med', 'Width_first','Rheobase', 'Rin', 'Cm',...,
    'AHP_amp','AHP_dur','AHP_area'}; 

data_temp = cell(1,numel(fields));
struct_temp = cell(1,size(dl,1));

for shi = 1:numel(fields)
    data_temp{shi} = xlsread(strcat(fp_import_file,'\',import_filename),shi);
end

for fi = 1:numel(fields)
    
    for sti = 1:size(dl,1)
        if data_temp{fi}(1,1) == 0 
            if size(data_temp{fi},2) < 5
                
                struct_temp{sti}.(fields{fi}) = data_temp{fi}(2:end,sti);

            else
                
                struct_temp{sti}.(fields{fi}) = data_temp{fi}(dl(sti,1):dl(sti,2),:);

            end
        else
            
               struct_temp{sti}.(fields{fi}) = data_temp{fi}(1:end,sti);

        end

    end
end

current_inj = (20:20:400)'; %in pA
shk3_DR_CNO = struct_temp{1};
shk3_CNO = struct_temp{2};
% WT = struct_temp{3};
% WT_saline = struct_temp{4};

%% save as mat data
if save_results == 1
    cd (fp_analyzed_data)

    save(filename,'shk3_CNO','shk3_DR_CNO','current_inj')
end
