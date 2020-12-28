function [timepoints] = plot_trace(filename,ch,n)
%plot_trace('filename') plots the data as a function of time
%   Take the nth sweep in a h5 data file and plot it as a function of time
%   The channel is defined by ch
%   Return the timestamp of each data point
%   Installation of WaveSurfer Required
if nargin<3           %default for the sweeps is a continous display of all sweeps
    n=0;
end

extracted_data = ws.loadDataFile(filename);
time_increment = 1/extracted_data.header.AcquisitionSampleRate;
fields = fieldnames(extracted_data);  %get the fieldnames of the data struct

if n == 0
    number_of_datapoints = size(extracted_data.(fields{n+2}).analogScans(:,ch),1); %get the number of data points of the sweep
else
    number_of_datapoints = size(extracted_data.(fields{n+1}).analogScans(:,ch),1);
end

endtime = time_increment*number_of_datapoints;
timepoints = (0:time_increment:endtime);
timepoints = timepoints(1:end-1);

total_filename = numel(filename); %number of characters in the filename

if contains(filename,'-')
    if n==0
    for i = 1:(length(fields)-1)
        figure('position',[56 200 1000 490]);
        plot(timepoints,extracted_data.(fields{i+1}).analogScans(:,ch),'k-','Linewidth',1);
        xlabel('Time (s)');
        xlim([0 endtime]);
        ylabel('mV');
        title(strcat(filename(1:total_filename-3),': ', fields(i+1)), 'Interpreter', 'none');
        
        waitforbuttonpress
    end
    else
    figure('position',[56 200 1000 490]);
    plot(timepoints,extracted_data.(fields{n+1}).analogScans(:,ch),'k-','Linewidth',1);
    xlabel('Time (s)');
    xlim([0 endtime]);
    ylabel('mV');
    title(strcat(filename(1:total_filename-3),': ', fields(n+1)), 'Interpreter', 'none');
    
    end
else
    figure('position',[56 200 1400 490]);
    plot(timepoints,extracted_data.(fields{n+2}).analogScans(:,ch),'k-','Linewidth',1);
    xlabel('Time (s)');
    xlim([0 endtime]);
    if endtime >=15
        ylabel('mV'); %for current clamp
    else
        ylabel('pA'); %for seal tests, y axis is in pA
    end
    
   title(strcat(filename(1:total_filename-3),': ', fields(n+2)), 'Interpreter', 'none');
end
