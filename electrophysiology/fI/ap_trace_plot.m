%% Select which cell and which trace to plot

% where to save the selected trace data
save_data_path = 'C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di\representative_traces\fI';

% name of the saved data
save_data_name = 'chronic_hm4di_fI_adult_ctrl_24_48.mat';

% analyzed fI results folder
fp_data = 'C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di\analyzed_fI_results';

% experimental conditions
cond = {'CNO','DR_24h','DR_48h'};

% Choose f_I analysis mat file
experiment = {'fI_200215', 'fI_200214','fI_201118'};

% Choose a cell in each experiment for plotting
cell_num = [8 1 3];

% Plotting mode (0 for non-normarlized traces, 1 for normalized traces) 
plot_mode = 0;

% current step
stim = 200;

%tint factors (btw 0 and 1, closer to 1 = more tint)
tint_factor = [0 0 0.5];

%tint and convert function
    %x-RGB value, y-tint factor
fun = @(x,y) (x+(255-x)*y)./255;

color1 = fun([0 0 0],tint_factor(1)); %black
color2 = fun([70 130 180],tint_factor(2)); %steel blue
color3 = fun([70 130 180],tint_factor(3)); %steel blue+50% tint

colorcode = cat(1, color1, color2, color3);
%% data readout
% Variable initializations
cell_data = cell(1,numel(experiment));
trace_data = cell(1,numel(experiment));
trace_id = NaN(1,numel(experiment));
rheo_id = NaN(1,numel(experiment));

for fi = 1:numel(experiment)
    curr_file = matfile(strcat(fp_data,'\',experiment{fi}, '.mat'));
    curr_file_data = curr_file.aDAT;
    curr_cell_id = curr_file.cell_id;
    curr_rheobase_ind = curr_file.rheobase_ind;
    
    cell_data{1,fi} = curr_file_data{1,cell_num(fi)};
    rheo_id(1,fi) = curr_rheobase_ind(fi,1);
    
    if plot_mode == 0
        trace_id(1,fi) = curr_cell_id{1,cell_num(fi)}(1,1)+stim/20-1;
    elseif plot_mode == 1
        trace_id(1,fi) = curr_cell_id{1,cell_num(fi)}(1,1)+rheo_id(1,fi)-1+stim/20;
    end
    
    trace_data{1,fi} = cell_data{1,fi}(:,trace_id(1,fi));
end

%% Plotting 

plot_range = 8000:22000;
figure('position',[56 200 1000 490])
hold on
for ci = 1:numel(cell_num)
    plot(trace_data{1,ci}(plot_range,1),'Color',colorcode(ci,1:3),'LineWidth',3)
end

%draw scale
plot([13500; 15500],[0; 0], '-k',[13500;13500],[0; 20], '-k', 'LineWidth',2)
text(13400, 5, '20 mV', 'HorizontalAlignment', 'right')
text(14000, -5, '0.2 s', 'HorizontalAlignment', 'center')

legend(cond)

if plot_mode == 0
    title(strcat(num2str(stim), 'pA injected'))
elseif plot_mode == 1
    title(strcat(num2str(stim),'pA over rheobase'))
end

hold off

%% Drawing current step

figure('position',[56 200 1000 490])
bl_before = zeros(1,1999);
current_step = repmat(stim/1000, 1, 10000);
bl_after = zeros(1,1999);

current_all = horzcat(bl_before,current_step,bl_after);

plot(current_all,'k', 'LineWidth',3)
ylim([-0.2 1])
    
%% save data
cd(save_data_path)

save(save_data_name, 'cond','experiment','cell_num','trace_data')
