%% groups for plotting
%name of groups
groups = {'adult_CNO_24_n_48','adult_DR_CNO','adult_DR_CNO_48h'};

%label for each group
g_label = {'CNO','DR+24h','DR+48h'};

%name of the variable 
parameter = 'Lat';

%name of the varibale for label
parlabel = 'Latency';

%unit of the variable
unit = 'ms';

%scale for unit conversion
scale = 1;

%exponent value of the Y axis
y_axis_expo = 0;

%lower and upper limit of the Y axis
y_limit = [-inf inf]; %inf indicates auto [-inf,inf]

%direction of the Y axis
y_dir_all = {'normal','reverse'};
y_dir = y_dir_all{1};

%plotting mode: 0 for single column variable, 1 for multiple column
%variable
plot_mode = 1;

%data selection mode: only applicable when plot_mode == 1
%1 for no normalization(start from the 1st current step)
%2 for normalization (start from the rheobase step)
select_mode = 1;

%stimulation amplitude: only applicable when plot_mode == 1
%indicates the stimulation amplitude (current injected, in pA)
stim = 300;

%offset of scatter points with respect to the error bar
dp_offset = 0.2;

%tint factors (btw 0 and 1, closer to 1 = more tint)
tint_factor_1 = [0 0 0.5]; %bar tint
tint_factor_2 = [0 0 0]; %point tint

%tint and convert function
    %x-RGB value, y-tint factor
fun = @(x,y) (x+(255-x)*y)./255;

%col1-3 store face color, col4-6 store edge color
%bar

col_pp1 = 1:3;
col_pp2 = 4:6;

col_bar = NaN(numel(groups),3);
% col_bar(1,col_pp1) = fun([255 203 164]); %deep peach
% col_bar(2,col_pp1) = fun([153 186 221]); %carolina blue
% col_bar(3,col_pp1) = fun([144 200 144]); %light green
col_bar(1,col_pp1) = fun([255 255 255],tint_factor_1(1)); %white
col_bar(2,col_pp1) = fun([70 130 180],tint_factor_1(2)); %salmon
col_bar(3,col_pp1) = fun([70 130 180],tint_factor_1(3)); %steel blue
% col_bar(3,col_pp1) = fun([250 128 114],tint_factor_1(3)); %salmon
% col_bar(4,col_pp1) = fun([70 130 180],tint_factor_1(4)); %steel blue

% col_bar(1,col_pp1) = fun([255 203 164]); %deep peach
% col_bar(2,col_pp1) = fun([153 186 221]); %carolina blue
%col_bar(3,1:3) = fun([255 203 164]); %deep peach
%col_bar(4,1:3) = fun([153 186 221]); %carolina blue

% col_bar(1,col_pp2) = fun([255 203 164]); %deep peach
% col_bar(2,col_pp2) = fun([153 186 221]); %carolina blue
% col_bar(3,col_pp2) = fun([144 200 144]); %light green
col_bar(1,col_pp2) = fun([0 0 0],tint_factor_1(1)); %black
col_bar(2,col_pp2) = fun([70 130 180],tint_factor_1(2)); %salmon
col_bar(3,col_pp2) = fun([70 130 180],tint_factor_1(3)); %steel blue
% col_bar(3,col_pp2) = fun([250 128 114],tint_factor_1(3)); %salmon
% col_bar(4,col_pp2) = fun([70 130 180],tint_factor_1(4)); %steel blue

% col_bar(1,col_pp2) = fun([255 203 164]); %deep peach
% col_bar(2,col_pp2) = fun([153 186 221]); %carolina blue
%col_bar(3,4:6) = fun([255 203 164]); %deep peach
%col_bar(4,4:6) = fun([153 186 221]); %carolina blue


%points
col_poi = NaN(numel(groups),3);
col_poi(1,col_pp1) = fun([0 0 0],tint_factor_2(1)); %black
col_poi(2,col_pp1) = fun([211 211 211],tint_factor_2(2)); %light gray
col_poi(3,col_pp1) = fun([105 105 105],tint_factor_2(3)); %dim gray
% col_poi(3,col_pp1) = fun([211 211 211],tint_factor_2(3)); %light gray
% col_poi(4,col_pp1) = fun([105 105 105],tint_factor_2(4)); %dim gray
% col_poi(1,1:3) = fun([211 211 211]); %light gray
% col_poi(2,1:3) = fun([105 105 105]); %dim gray

col_poi(1,col_pp2) = fun([0 0 0],tint_factor_2(1)); %black
col_poi(2,col_pp2) = fun([211 211 211],tint_factor_2(2)); %light gray
col_poi(3,col_pp2) = fun([105 105 105],tint_factor_2(3)); %dim gray
% col_poi(3,col_pp2) = fun([211 211 211],tint_factor_2(3)); %light gray
% col_poi(4,col_pp2) = fun([105 105 105],tint_factor_2(4)); %dim gray
% col_poi(1,4:6) = fun([211 211 211]); %light gray
% col_poi(2,4:6) = fun([105 105 105]); %dim gray

%transparency of the face color
falpha = [1 1 1 1];
%transparency of the edge color
ealpha = [1 1 1 1];


%% prepare data for plotting
data = cell(1,numel(groups));
ave_data = NaN(1,numel(groups));
std_data = NaN(1,numel(groups));
sem_data = NaN(1,numel(groups));
er = cell(1,numel(groups));
X_rep = cell(1,numel(groups));
X_swarm = cell(1,numel(groups));
Y_swarm = cell(1,numel(groups));

% positions on the X axis
X = (1:numel(groups));

% data readout and stats calculation
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
    
    X_rep{1,gi} = repmat(X(1,gi),1,numel(data{1,gi}));
    [X_swarm{1,gi}, Y_swarm{1,gi}] = swarmplot(X_rep{1,gi}+dp_offset,data{1,gi},0.1);
end

%% plotting

if numel(groups) <=3
    figure()
else
    figure('Position',[500 100 700 500]);
end

for pi = 1:numel(groups)
    
    bar(X(1,pi), ave_data(1,pi),... 
        'EdgeColor', col_bar(pi,4:6), ...
        'EdgeAlpha', ealpha(pi),...
        'FaceColor', col_bar(pi,1:3), ...
        'FaceAlpha', falpha(pi),...
        'LineWidth',2, 'BarWidth', 0.4)
    hold on

    er{1,pi} = errorbar(X(1,pi), ave_data(1,pi), sem_data(1,pi));
    er{1,pi}.Color = 'k';
    er{1,pi}.LineWidth = 2;
    er{1,pi}.CapSize = 12;
    
    scatter(X_swarm{1,pi},Y_swarm{1,pi},'MarkerEdgeColor',col_poi(pi,4:6),...
        'MarkerFaceColor',col_poi(pi,1:3),'LineWidth',1,'SizeData',20)
end

%aesthetics of the plot
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 2;
ax.YLabel.String = strcat(parlabel,' (',unit,')');
ax.YAxis.Exponent = y_axis_expo;
ax.YDir = y_dir;

if sum(isinf(y_limit)) ~= 2
    ax.YLim = y_limit;
end
%ax.YTick = [0 50 100];

ax.XTick = X;
%ax.XTick = [];
% ax.XTickLabel = [];
ax.XTickLabel = g_label;
ax.TickLabelInterpreter = 'none';
ax.XAxisLocation = 'origin';
%ax.XTickLabel.FontSize = 14;
%ax.XTickLabelRotation = 60;
%legend('hM4Di+CNO','Scaled','ctrl+CNO','FontSize',12,'Location','southeast')

hold off
box off
