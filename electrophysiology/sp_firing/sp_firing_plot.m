%%%%This script takes a number of h5 files recording spontaneous firing 
%%%%of a neuron and concacate them into a single plot. The average firing
%%%%rate of each file is also calculated.

%%%%% h5 files should have more than 30 s of recording
%%%%% in the current setting, 5 control file followed by 5 CNO file

%%%% whether to show plot
figure_on = 1;

%%%% whether to save to file
save_on = 0;

all_files = dir; %files in the current folder
file_num = numel(all_files)-2; %all data files
sprate = 5000; %rate after resampling

filename = cell(1,file_num);

cell_id = cell(1,file_num);
% for cjj = 1:file_num
%     cell_id{1,cjj} = zeros(file_num,1); %first column cell ID, second column trace ID
% end

aDAT = cell(1,file_num);
for cii = 1:file_num
    aDAT{1,cii} = zeros(30*sprate,file_num); %pre-allocate each cell in the aDAT cell array
end

%Get file names of the specified cell
for f = 1:file_num
    filename{f} = all_files(f+2).name;
    num_filename = numel(filename{f}); %number of characters in the file name
    
    
    cellID = str2double(filename{f}(5));
    trace_id = str2double(filename{f}(num_filename-6 : num_filename-3));
    
    
    cell_id{1,cellID}(trace_id,1) = trace_id;
    
    
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
    dtstruct = ws.loadDataFile(filename{f});
    prefix = strcat('dtstruct.sweep_',sname,'.analogScans');
    vals = eval(prefix);
    vals = vals(:,1);
    
    %resample @1:2 (sample rate = 5000 Hz)
    resampled(1:numel(resample(vals,1,2)),1) = resample(vals,1,2);
    
    %Cut data at 30 seconds
    %In aDAT, each cell represents a neuron; within each cell,each column
    %is a single file of 30 seconds
    aDAT{1,cellID}(1:30*sprate,trace_id) = resampled(1:30*sprate);
end

for gi = 1:cellID
    cell_id{1,gi} = nonzeros(cell_id{1,gi}); %remove zeros in the cell_id array
end

%%%%Calculating firing rate for each file (same cell)
num_of_spikes = zeros(10,cellID); %number of spikes
FR = zeros(10,cellID); %store firing rates
Vthresh = zeros(10,cellID); %store spiking thresholds
aveFR = zeros(cellID,2); %store average baseline and experimental firing rates

V_thresh = -20;

for ci = 1:cellID
    
    if isempty(cell_id{1,ci}) == 1
        continue
    else
       fi_start = cell_id{1,ci}(1,1);
    end
        
        for fi = fi_start:fi_start+10-1
            if fi > size(aDAT{1,ci},2)
                continue
            else
                filedata = aDAT{1,ci}(:,fi);
            end
            
            Vthresh(fi,ci) = get_spike_threshold(filedata,20);

            for ti = 2:30*5000
                if filedata(ti)>=V_thresh && filedata(ti-1)<V_thresh
                    num_of_spikes(fi-fi_start+1,ci)= num_of_spikes(fi-fi_start+1,ci)+1;
                end
            end

            FR(fi-fi_start+1,ci) = num_of_spikes(fi-fi_start+1,ci)/30; %firing rates in Hz
        end

        aveFR(ci,1) = mean(FR(1:5,ci)); %firing rate of baseline
        aveFR(ci,2) = mean(FR(6:10,ci)); %firing rate after CNO onset
end


%%%%Putting all file together and plotting
concatDAT = zeros(150000*10,cellID); 
timestamp = zeros(150000*10,cellID);
%%    
for cj = 1:1%cellID
    
    if isempty(cell_id{1,cj}) == 1
        continue
    else
       i_start = cell_id{1,cj}(1,1);
    end
    
    for i = i_start:i_start+10-1 
       %each cell should have 10 files: 5 for baseline and 5 for CNO
       %change the number of loops if # of files per cell/condition is
       %different
       if i > size(aDAT{1,cj},2)
                continue
       else

            concatDAT((150000*(i-i_start)+2):150000*(i-i_start+1),cj)= ...
               aDAT{1,cj}(2:end,i);
       end
    end
    
    for ni = 1:150000*10
        if concatDAT(ni,cj) == 0
            concatDAT(ni,cj) = NaN;
        end
    end


%%%% Generating timestamps    
    for j = 1:size(concatDAT(:,cj),1)
        timestamp(j,cj) = j/5000;
    end
    

    if figure_on == 1
        figure('position',[56 200 1400 300])
        plot(timestamp(:,cj),concatDAT(:,cj),'k-');
        xlabel('Time (s)')
        ylabel('mV')
        title(strcat('Cell',num2str(cj)))
        set(gca,'ylim',[-80 40]);
        box off
    end
  
end

%%%%Save to file
if save_on == 1
    current_dir_name = cd; 
    total_dir_name_number = numel(current_dir_name);
    display_name = strcat(current_dir_name(total_dir_name_number-5: end),'.mat');
    save(display_name,'num_of_spikes','FR','aveFR','Vthresh')
    
else
    disp('Data not saved!')
end