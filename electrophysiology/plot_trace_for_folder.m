%%%%Plot traces for a folder

all_files = dir;
file_num = numel(all_files)-2;

filename = cell(1,file_num);
filename_bycell = cell(1,file_num);
cell_id = cell(1,file_num);
cell_num = NaN(1,file_num);
exDAT = cell(1,file_num);
timepoints = cell(1,file_num);



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
    
    if trace_id < 10
                sname = strcat('000',num2str(trace_id));
            elseif trace_id >= 10 && trace_id < 100
                sname = strcat('00',num2str(trace_id));
            elseif trace_id >= 100 && trace_id < 1000
                sname = strcat('0',num2str(trace_id));
            else
                sname = num2str(trace_id);
    end
    
    dtstruct = ws.loadDataFile(filename{f});
    prefix = strcat('dtstruct.sweep_',sname,'.analogScans');
    vals = eval(prefix);
    exDAT{1,cellID}(1:numel(vals(:,1)),trace_id) = vals(1:numel(vals(:,1)),1); 
    
    num_of_datapoints = size(vals,1);
    
    time_increment = 1/dtstruct.header.AcquisitionSampleRate;
    endtime = time_increment*num_of_datapoints;
    current_timepoints = (0:time_increment:endtime);
    current_timepoints = current_timepoints(1:end-1);
    timepoints{1,cellID}(1:numel(current_timepoints),trace_id) = current_timepoints;
    
end

for ii = 1:max(cell_num)
    if isempty(cell_id{1,ii}) == 1
        continue
    end
    
    for jj = 1:cell_id{1,ii}(end)
        
    temp = nonzeros(exDAT{1,ii}(:,jj));
    exDAT{1,ii}(:,jj) = NaN;
    exDAT{1,ii}(1:numel(temp),jj) = temp;
    end
end

for i = 1:max(cell_num)
    if isempty(cell_id{1,i}) == 1
        continue
    end
    
    figure('position',[56 200 1000 490]);
    sgtitle(strcat('Cell',num2str(i)))
    for j = 1:size(cell_id{i},1)
        if cell_id{1,i}(j,1) == 0
            continue
        end
        
        
        subplot(5,4,j);
        plot(timepoints{1,i}(:,j),exDAT{1,i}(:,j),'k')
        xlabel('Time (s)');
        title(filename_bycell{1,i}(j),'Interpreter','none');
    end
end
    