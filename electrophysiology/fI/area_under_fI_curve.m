%% import mat files, areas stored by variables, within each variable, each column indicates an experimental condition.

%location where the mat file will be saved
fp_analyzed_data = 'C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di\area_under_curve_fI';

%name of the saved file
filename = 'Area_Under_Curve_adult_fI.mat';

%save results
save_file = 1;

% names of experimental conditions within the analyzed fI mat file
condn = {'WT_CNO','DR_CNO'};

adult_CNO = struct;
adult_DR_CNO = struct;

%field names
fields = {'MFR', 'IFR', 'mean_IFR', 'Threshold', 'ADP_ind','Lat',...
    'Udratio','Width'}; 

%field names
% field1 = 'MFR'; 
% field2 = 'IFR';
% field3 = 'mean_IFR';
% field4 = 'Threshold';
% field5 = 'ADP_ind';
% field6 = 'Lat';
% field7 = 'Udratio';
% field8 = 'Width';

%% Parameter initialization
cell_num = NaN(numel(condn),1);

for gi = 1:numel(condn)
    cell_num(gi,1) = size(eval(strcat(condn{gi},'.',field1)),2);
end

data_temp = NaN(max(cell_num),numel(condn));
A_temp = NaN(max(cell_num),numel(condn));

AreaUnderCurve = struct(field1, data_temp,...
    field2, data_temp,...
    field3, data_temp,...
    field4, data_temp,...
    field5, data_temp,...
    field6, data_temp,...
    field7, data_temp,...
    field8, data_temp);

fn = fieldnames(AreaUnderCurve);

%% calculate areas under each individual curve and save in the struct

%each gi corresponds to an experimental condition

for gi = 1:numel(condn) %per condition
   for cii = 1:max(cell_num) %per cell
       for fii = 1:numel(fn) %per field
           current_data = eval(strcat(condn{gi},'.',fn{fii}));
           if isnan(current_data(1:20,cii))
               A_temp(cii,gi) = NaN;
           else
               A_temp(cii,gi) = areaundercurve(current_inj,current_data(1:20,cii));
           end
           
           AreaUnderCurve.(fn{fii})(cii,gi) = A_temp(cii,gi);
       end
   end
end


%% save to file
if save_file == 1
    %location where the mat file will be saved
    fp_analyzed_data = 'C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di';
    filename = {'Areas_Under_Curve_fI_adult.mat'};

    cd(fp_analyzed_data)
    save(filename{1},'AreaUnderCurve')
end