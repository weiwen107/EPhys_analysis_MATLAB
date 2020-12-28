%% these need to be changed for each new data run

%Name of file you would like to have data saved to
save_file_name = 'SPfiring_200809.mat';

%loop through all the folders
experiment = '200809';

%Filepath where you keep data folders
fp_data = 'C:\Users\schum\Google_Drive\Lab\Data\acute_hm4di\sp_firing\';

%location to save the analyzed mini data
fp_analyzed_data = ...,
    'C:\Users\schum\Google_Drive\Lab\Data_analysis\acute_DREADDs\acute_hm4di\analyzed_sp_data\';

%% 
% whether to save to file
save_results = 1;

%sampling rate for analysis
sp_rate = 5000;

%AP detection threshold
V_thresh = -20;
dV_thresh = 20;

current_folder = strcat(fp_data,experiment);
cd(current_folder)

all_files = dir; %files in the current folder
file_num = numel(all_files)-2; %all data files
filename = cell(1,file_num);
cell_id = cell(1,file_num);
cell_num = NaN(1,file_num);
filename_bycell = cell(1,file_num);

aDAT = cell(1,file_num);
num_of_spikes = cell(1,file_num); %number of spikes
FR = cell(1,file_num); %store firing rates
Vthresh = cell(1,file_num); %store spiking thresholds


%Create cell id index
for f = 1:file_num
    filename{f} = all_files(f+2).name;
    num_filename = numel(filename{f}); %number of characters in the file name
    
    
    if isnan(str2double(filename{f}(6))) == 1
        cellID = str2double(filename{f}(5));
    else
        cellID = str2double(filename{f}(5:6));
    end
    
    cell_num(f) = cellID;
    trace_id = str2double(filename{f}(num_filename-6 : num_filename-3));


    cell_id{1,cellID}(trace_id,1) = trace_id;

    filename_bycell{1,cellID}{trace_id} = filename{f};
end

for gi = 1:max(cell_num)
    cell_id{1,gi} = nonzeros(cell_id{1,gi});
    if isempty(cell_id{1,gi}) == 1
        continue
    else
        for ri = 1:cell_id{1,gi}(end,1)
            if ismember(ri,cell_id{1,gi}) == 0
                aDAT{1,gi}(:,ri) = NaN;
            end
        end
    end
end

%% Data readout and FR calculation

for ci = 1:max(cell_num)
    if isempty(cell_id{1,ci}) == 1
            continue
    else
        for ti = 1:size(cell_id{1,ci},1)
                trace_id = cell_id{1,ci}(ti,1);
    
            if trace_id < 10
                        sname = strcat('000',num2str(trace_id));
                    elseif trace_id >= 10 && trace_id < 100
                        sname = strcat('00',num2str(trace_id));
                    elseif trace_id >= 100 && trace_id < 1000
                        sname = strcat('0',num2str(trace_id));
                    else
                        sname = num2str(trace_id);
            end
    
    
            %datafile readout
            dtstruct = ws.loadDataFile(filename_bycell{1,ci}{trace_id});
            prefix = strcat('dtstruct.sweep_',sname,'.analogScans');
            vals = eval(prefix);
            vals = vals(:,1);

            %resample @1:2 (sample rate = 5000 Hz)
            resampled(1:numel(resample(vals,1,2)),1) = resample(vals,1,2);
            aDAT{1,ci}(1:numel(resampled(:,1)),trace_id) = resampled;
            
            clear resampled

            %Calculating firing rate for each file
            filedata = aDAT{1,ci}(:,trace_id);
            time = numel(filedata)/sp_rate; % in seconds

            sp_counter = 0;

            for tii = 1:numel(filedata)-1
                if filedata(tii)>=V_thresh && filedata(tii+1)<V_thresh
                    sp_counter = sp_counter+1;
                end
            end
            
            num_of_spikes{1,ci}(ti,1) = sp_counter;
            FR{1,ci}(ti,1) = num_of_spikes{1,ci}(ti,1)/time; %firing rates in Hz
            
            if sp_counter > 0
                Vthresh{1,ci}(ti,1) = get_spike_threshold(filedata,dV_thresh);
            else
                Vthresh{1,ci}(ti,1) = NaN;
            end
        end
    end
end

%% save results
if save_results == 1
    cd(fp_analyzed_data)
    save(save_file_name,'cell_id', 'Vthresh', 'num_of_spikes','FR')
end