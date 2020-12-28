sample1 = data{1};
sample2 = data{2};
sample3 = data{3};
% sample4 = shk3_WT_CNO_DW.decaytau;
% sample5 = shk3_WT_DR_CNO_DW.decaytau;
% sample5 = shk3_CNO.amp;
% sample6 = shk3_DR_CNO.amp;

%choose the type of test: 0 for onw-way ANOVA, 1 for kruskal-wallis
p = 1;

%choose the type of post-hoc test
ph = {'tukey-kramer','bonferroni','dunn-sidak','lsd','scheffe'};

%prepare data set for the test
y = [sample1' sample2' sample3'];% sample4' sample5'% sample6'];
g1 = repmat('A',1,numel(sample1));
g2 = repmat('B',1,numel(sample2));
g3 = repmat('C',1,numel(sample3));
% g4 = repmat('D',1,numel(sample4));
% g5 = repmat('E',1,numel(sample5));
% g6 = repmat('F',1,numel(sample6));
g = [g1 g2 g3];% g4 g5];
g_cell = num2cell(g);

if p == 0

    %one-way ANOVA
    [p, tbl, stats] = anova1(y,g_cell);
    [c, m, h, hms] = multcompare(stats, 'alpha', .05, 'ctype', ph{1});
    
else
    
    %Kruskal-wallis
    [p, tbl, stats] = kruskalwallis(y,g_cell);
    [c, m, h, hms] = multcompare(stats, 'alpha', .05, 'ctype', ph{1});
    
end