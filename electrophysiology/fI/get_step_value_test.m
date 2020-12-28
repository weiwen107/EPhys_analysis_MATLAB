% input data (only applies to processed data, not raw data)
data = WT_CNO.ADP_ind;

% mode=1: no normalization (starting from the 1st current step)
% mode=2: normalized to the rheobase (starting from rheobase step)
mode = 2;

% Rheobase
rheobase = WT_CNO.Rheobase;

%stimulation amplitude (current injected, in pA)
stim = 200; 

cell_num = size(data,2); % how many cells
trace_num = size(data,1); % how many current steps

% output data
vals = NaN(cell_num,1);

rheo_ind = NaN(1,cell_num);
trace_ind = NaN(1,cell_num);

for ci = 1:cell_num
    

    rheo_ind(1,ci) = rheobase(ci,1)/20;

    if mode == 1
        if stim == 0
            warning('Trace 0 does not exist!')
        else
            trace_ind(1,ci) = stim/20;
        end
    elseif mode == 2
        trace_ind(1,ci) = rheo_ind(1,ci) + stim/20;
    else
        warning('Unknown mode')
    end

    if isnan(trace_ind(1,ci)) || trace_ind(1,ci) > 20
        vals(ci,1) = NaN;
        warning(strcat('Cell', num2str(ci),' Trace', num2str(trace_ind(1,ci)),' does not exist!'))
    else

        vals(ci,1) = data(trace_ind(1,ci),ci);

    end
    
end