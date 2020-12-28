function [areasum] = areaundercurve(x_data,y_data,plot_on)
% Given a dataset, calculate the area under the curve by dividing the
% points into trapezoids

if nargin < 3
    plot_on = 0;
end

% x_data = current_inj;
% y_data = DR_CNO.MFR(:,1);
trapi = NaN(numel(y_data)-1,1);

if numel(x_data) ~= numel(y_data)
    print('X and Y should have same size')
else
    for ti = 1:size(trapi,1)
        if isnan(y_data(ti)) || isnan(x_data(ti))
            continue
        else
            trapi(ti) = (y_data(ti)+y_data(ti+1))*(x_data(ti+1)-x_data(ti))*0.5;
        end
    end
    areasum = nansum(trapi);
    
   if plot_on
        figure;
        plot(x_data,y_data,'o-')
        hold on
        for tii = 1:size(y_data,1)
            if isnan(y_data(tii))
                continue
            else
                plot([x_data(tii) x_data(tii)],[0 y_data(tii)],'r-')
            end
        end
    end
end
end