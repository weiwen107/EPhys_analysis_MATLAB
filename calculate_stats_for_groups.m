%% groups for stats calculation
%name of groups
groups = {'shk3_CNO','shk3_DR_CNO'};

%name of the variable 
parameter = 'Rin';

%scale for unit conversion
scale = 1;

%plotting mode: 0 for single column variable, 1 for multiple column
%variable
plot_mode = 0;

%data selection mode: only applicable when plot_mode == 1
%1 for no normalization(start from the 1st current step)
%2 for normalization (start from the rheobase step)
select_mode = 1;

%stimulation amplitude: only applicable when plot_mode == 1
%indicates the stimulation amplitude (current injected, in pA)
stim = 300;

%% Stats calculation

data = cell(1,numel(groups));
ave_data = NaN(1,numel(groups));
std_data = NaN(1,numel(groups));
sem_data = NaN(1,numel(groups));

for gi = 1:numel(groups)
    if plot_mode == 0
    
        data{1,gi} = eval(strcat(groups{gi},'.',parameter));
        data{1,gi} = data{1,gi}.*scale;
    elseif plot_mode == 1
       
        data_temp = eval(strcat(groups{gi},'.',parameter));
        rheo = eval(strcat(groups{gi},'.','Rheobase'));
        data{1,gi} = get_step_values(data_temp, rheo, select_mode, stim);
    end
    
    
    ave_data(1,gi) = nanmean(data{1,gi});
    std_data(1,gi) = nanstd(data{1,gi});
    sem_data(1,gi) = nansem(data{1,gi});
    
end
