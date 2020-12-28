function [V_thresh] = get_spike_threshold(data, dV_thresh, plot_on)

if nargin < 3
    plot_on = 0;
end

V = data;
dV = diff(V);

% get dV/dt (convert mV/s to V/s)
dV_sec = dV.*10;
% find when dV/dt crosses 20 V/s threshold
fourth_sample = find(dV_sec > dV_thresh,4,'first'); %used to weed out intial dV spike (first 4)
if numel(fourth_sample) < 4
    V_thresh = NaN;
    V_xi = NaN;
else
next_sample = fourth_sample(4);
prev_sample = find(dV_sec(1:next_sample) < dV_thresh,1,'last');

% linear interpolation to find exact sample point at which the cross
% happened
dV_xi = dV_thresh;
dV_x = [dV_sec(prev_sample) dV_sec(next_sample)]'; %dV/dt
dV_y = [prev_sample next_sample]'; % index of dV/dt


    try
        dV_yi = interp1q(dV_x,dV_y,dV_xi);
    catch
        keyboard;
    end
    
    % linear interpolation on V to find Vm at which dV/dt = 20 V/s
    V_x = dV_y + 1;
    V_y = V(V_x);
    V_xi = dV_yi + 1;
    V_yi = interp1q(V_x,V_y,V_xi);
    
    % this is our threshold
    V_thresh = V_yi;
end
%%
if plot_on
    figure(); hold on;
    plot(V,'linewidth',2);
    plot(dV_sec./10,'linewidth',1.5);
    plot(V_xi,V_thresh,'o','markersize',10,'markerfacecolor','k');
    set(gca,'xlim',[8000 15000]);
end



