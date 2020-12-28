%data readout
dtstruct = ws.loadDataFile('cell2_0004.h5');
vals = dtstruct.sweep_0004.analogScans;
vals = vals(:,1);

%resample @ 1:2
resampled(1:numel(resample(vals,1,2)),1) = resample(vals,1,2);

%choose the range of data for plotting
% f = figure('position',[119 171 1210 611]);
% hold on
% axis([1 (numel(resampled)+100) (nanmean(resampled)-20) (nanmean(resampled)+100)])

% plot(resampled)
% 
% [X,Y] = ginput;
% t_start = int32(X(1));
% t_end = int32(X(2)); 
% hold off
% close(f)

%convert x coordinates to minutes
timestamp = (1:numel(resampled))*0.2/1000/60;

figure('position',[119 171 1210 611])
plot(timestamp,resampled,'k');
ylim([-70 40])
xlabel('Time (min)')
title(dtstruct.header.DataFileBaseName)