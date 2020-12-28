%% groups and parameters for testing (no. of groups >=3)
%name of groups
groups = {'adult_CNO_24_n_48','adult_DR_CNO','adult_DR_CNO_48h'};

%name of the variable 
parameter = 'Rin';

%testing mode: 0 for single column variable, 1 for multiple column
%variable
test_mode = 0;

%data selection mode: only applicable when plot_mode == 1
%1 for no normalization(start from the 1st current step)
%2 for normalization (start from the rheobase step)
select_mode = 1;

%stimulation amplitude: only applicable when plot_mode == 1
%indicates the stimulation amplitude (current injected, in pA)
stim = 300;

%choose the type of post-hoc correction
ph = {'tukey-kramer','bonferroni','dunn-sidak','lsd','scheffe'};
ph_curr = ph{1};

%confidence level
cl = 0.05;

%% data readout and preparation for multi-group tests

sample = cell(1,numel(groups));
counter = NaN(1,numel(groups));
counter_sum = 0;
g_alph = cell(1,numel(groups));

for gi = 1:numel(groups)
    
    if test_mode == 0
        sample{1,gi} = eval(strcat(groups{gi},'.',parameter));
    elseif test_mode ==1
        data_temp = eval(strcat(groups{gi},'.',parameter));
        rheo = eval(strcat(groups{gi},'.','Rheobase'));
        sample{1,gi} = get_step_values(data_temp, rheo, select_mode, stim);
    end
    
    counter(1,gi) = size(sample{1,gi},1);
    counter_sum = counter_sum + counter(1,gi);
    
    g_alph{1,gi} = repmat(char(gi+64),1,numel(sample{1,gi}));
end

y = NaN(1, counter_sum);
g = [];

start_pos = 1; 
end_pos = counter(1,1);

for yi = 1:numel(groups)
    y(1,start_pos : end_pos) = (sample{1,yi})';
    
    g = cat(2,g,g_alph{1,yi});
    
    if (yi + 1) <= numel(groups) 
        start_pos = start_pos + counter(1,yi);
        end_pos = end_pos + counter(1,yi+1);
    else
        continue
    end
    
end

g_cell = num2cell(g);

%% conducting test
% datasets undergo Shapiro-Wilk test first to see if they follow normal
% distribution. If all datasets show normal distribution, one-way ANOVA
% test will be used. If one dataset does not follow normal distribution,
% then the non-parametric Kruskal-wallis test will be used.

h_sw = NaN(1,numel(groups));
p_sw = NaN(1,numel(groups));
for hi = 1:numel(groups)
    [h_sw(1,hi), p_sw(1,hi), ~] = swtest(sample{1,hi},cl);
end

if sum(h_sw) == 0

    %one-way ANOVA
    [p, tbl, stats] = anova1(y,g_cell);
    [c, m, h, hms] = multcompare(stats, 'alpha', cl, 'ctype', ph_curr);
    
else
    
    %Kruskal-wallis
    [p, tbl, stats] = kruskalwallis(y,g_cell);
    [c, m, h, hms] = multcompare(stats, 'alpha', cl, 'ctype', ph_curr);
    
end