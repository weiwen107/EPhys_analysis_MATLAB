function [rise_avg,decay_avg] = get_avg_mini_kinetics(rise,decay)
%This function takes two cell arrays produced by the MINIANALYSIS script-
%RISE and DECAY, which store the rise time and decay time, respectively, of each trace of a specific
%cell.

rise_avg = NaN(20,numel(rise{1,1}));
decay_avg = NaN(20,numel(decay{1,1}));

if ~iscell(rise) || ~iscell(decay)
    warning('Inputs should be cell arrays.')
    
else
    
    for ci = 1:numel(rise{1,1})
        
        for ti = 1:size(rise{1,1}{1,ci},2)
            temp_1 = nonzeros(rise{1,1}{1,ci}(:,ti));
            rise{1,1}{1,ci}(:,ti) = NaN;
            rise{1,1}{1,ci}(1:numel(temp_1),ti) = temp_1;

            rise_avg(ti,ci) = nanmean(rise{1,1}{1,ci}(:,ti));
            
        end
    end
    
    for cj = 1:numel(decay{1,1})
        for tj = 1:size(decay{1,1}{1,cj},2)
            
            temp_2 = nonzeros(decay{1,1}{1,cj}(:,tj));
            decay{1,1}{1,cj}(:,tj) = NaN;
            decay{1,1}{1,cj}(1:numel(temp_2),tj) = temp_2;
            
            decay_avg(tj,cj) = nanmean(decay{1,1}{1,cj}(:,tj));
            
        end
    end

end

