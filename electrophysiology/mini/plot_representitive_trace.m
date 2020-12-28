%% Select which cell and which trace to plot

% where to save the selected trace data
save_data_path = 'C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di\representative_traces\mini';

% name of the saved data
save_data_name = 'chronic_hm4di_adult_DR_48CNO.mat';

% data folder
fp_data = 'C:\Users\schum\Google_Drive\Lab\Data\chronic_hm4di\mini\';

experiment = '201028';

trace_file = 'cell1_0002.h5';

%% data readout and plotting
extracted_data_file = ws.loadDataFile(strcat(fp_data,experiment,'\',trace_file));
fields = fieldnames(extracted_data_file);  %get the fieldnames of the data struct
sprate = extracted_data_file.header.AcquisitionSampleRate;
selected_trace_data = extracted_data_file.(fields{2}).analogScans(:,1);

time_increment = 1/sprate;
number_of_datapoints = size(selected_trace_data,1); %get the number of data points of the sweep

%pick the interval you want to show (in seconds)
start_time = 22;
startpoint = start_time * sprate;
end_time = 23;
endpoint = end_time * sprate;
timepoints = (start_time:time_increment:end_time);

plot_data = detrend(selected_trace_data(startpoint:endpoint,1));

figure('position',[56 200 700 490]);
plot(timepoints,plot_data,'k-','Linewidth',1);
hold on
xlabel('Time (s)');
xlim([start_time end_time]);
ylim([-100 50])
    
%draw scale
plot([start_time+0.8; start_time+0.9],[-80; -80], '-k',[start_time+0.8;start_time+0.8],[-80; -70], '-k', 'LineWidth',2)
text(start_time+0.79, -76, '10 pA', 'HorizontalAlignment', 'right')
text(start_time+0.85, -84, '0.1 s', 'HorizontalAlignment', 'center')
%set(gca, 'Visible', 'off')


%% save data
cd(save_data_path)
save(save_data_name, 'experiment','trace_file','start_time','end_time')