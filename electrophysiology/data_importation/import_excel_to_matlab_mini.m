%% import from excel
%location of the excel file
fp_import_file = 'C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di\matlab_import_file';
import_filename = 'matlab_import_mini_24.xlsx';

%location where the mat file will be saved
fp_analyzed_data = 'C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di\mini_data_by_groups';

%risetime cutoff
rise = 4;

%name of the mat file
filename = {strcat('chronic_DREADDs_48h_adult_rise',num2str(rise),'.mat')};

%save results
save_results = 1;

%experimental conditions,each corresponds to a sheet
cond = {strcat('adult_DR_CNO_48h_rise',num2str(rise)),strcat('adult_CNO_48h_rise',num2str(rise))};

%field names
fields = {'amp', 'frq', 'risetime', 'decaytau', 'charge',...
    'Vm', 'Ra', 'Cm', 'Rin'};

data_temp = cell(1, numel(cond));
struct_all =cell(1, numel(cond));

for shi = 1:numel(cond)
    data_temp{shi} = xlsread(strcat(fp_import_file,'\',import_filename),cond{shi});
end

for cdi = 1:numel(cond)
    
    for fi = 1:numel(fields)
        struct_all{cdi}.(fields{fi}) = data_temp{cdi}(:,fi);
    end
    
end

%% save data to corresponding data structure (by experimental condition)
adult_DR_CNO_48h = struct_all{1};
adult_CNO_48h = struct_all{2};
%DR_GFP_CNO = struct_all{3};

if save_results == 1
    cd(fp_analyzed_data)
    save(filename{1},'adult_DR_CNO_48h','adult_CNO_48h')
end