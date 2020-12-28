%%%% save data by properties rather than experimental conditions
% each property (e.g. amp) is stored as a table which includes all
% experimental conditions (in each column)

%% import from excel
%location of the excel file (also the place where mat file will be saved)
fp_analyzed_data = 'C:\Users\schum\Google_Drive\Lab\Data_analysis\Naspm';
import_filename = 'matlab_import_mini.xlsx';

%name of the mat file
filename = {'naspm_shk3_mini_by_property.mat'};

%save results
save_results = 1;

%experimental conditions (field names in this case)
cond = {'before','washin','after'}; %each corresponds to a sheet

%property names
properties = {'amp', 'frq', 'risetime', 'decaytau', 'charge',...
    'Vm', 'Ra', 'Cm', 'Rin'};

data_temp = cell(1, numel(cond));
cell_all =cell(1, numel(properties));

for shi = 1:numel(cond)
    data_temp{shi} = readtable(strcat(fp_analyzed_data,'\',import_filename),...,
        'ReadVariableNames',0,'FileType','spreadsheet','TreatAsEmpty','N/A','Sheet',cond{shi});
end

for cdi = 1:numel(properties)
    
    for fi = 1:numel(cond)
        cell_all{cdi}(:,fi) = data_temp{fi}(:,cdi);
    end
    
    cell_all{cdi}.Properties.VariableNames = cond;
end

%save data to corresponding data structure (by experimental condition)
amp = cell_all{1};
frq = cell_all{2};
risetime = cell_all{3};
decaytau = cell_all{4};
charge = cell_all{5};
Vm = cell_all{6};
Ra = cell_all{7};
Cm = cell_all{8};
Rin = cell_all{9};

%% save results
if save_results == 1
    cd(fp_analyzed_data)
    save(filename{1},'amp', 'frq', 'risetime', 'decaytau', 'charge',...
    'Vm', 'Ra', 'Cm', 'Rin')
end