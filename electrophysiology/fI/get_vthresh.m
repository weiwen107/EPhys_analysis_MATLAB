function [dV_sec, V_th, f_trough, sp_num, udratio] = get_vthresh(data, plot_on)
%%% This function calculates the threshold for AP by finding the point
%%% where dV/dt equals to 5% of the average maximum dV/dt from all APs
%%% V_th stores the threshold for each AP
%%% sp_index stores the indice of thresholds (first col) and fast trough
%%% (second col), indicating the starts and ends of each AP

if nargin < 2
    plot_on = 0;
end

V = data;
dV = diff(V);
dV_sec = smooth(dV.*10); %smoothed dV/dt (V/s)

sp_num = 0;
sp_peak = NaN(20,1); %store positions of AP peaks
sp_index = NaN(20,2); %start and end points of each spike
V_th = NaN(20,2); %threshold for each AP
f_trough = NaN(20,2); %fast trough for each AP
max_dv = NaN(20,1);
udratio = NaN(20,1); %upstroke/downstroke ratio

% count spikes
for di = 10000:20000
    if dV_sec(di) > 20 && dV_sec(di+1) < 20 && V(di) > -30
        sp_num = sp_num + 1;
    
        %if dV_sec(di)>0 && dV_sec(di+1)<15 %for each identified spike, find its peak and interval indices
            sp_peak(sp_num) = di;
            
            sp_index(sp_num,1) = find(dV_sec(10000:di) <= 5,1,'last')+10000;
            %sp_index(sp_num,2) = find(dV_sec(di+10:di+100)>=0.01*min(dV_sec(di+1:di+100)),1,'first')+di+1;
            %if dV/dt doesn't decrease to 1% of downstroke within 10 ms,
            %pick the maximum value within 5 ms after the peak.
            if isempty(find(dV_sec(di+10:di+100)>=0.01*min(dV_sec(di+1:di+100)),1,'first'))
                %[ft_val,ft_ind] = max(dV_sec(di+10:di+100));
                [ft_val,ft_ind] = min(V(di+10:di+80));
                sp_index(sp_num,2) = ft_ind+di+1;
            else
                sp_index(sp_num,2) = find(dV_sec(di+10:di+100)>=0.01*min(dV_sec(di+1:di+100)),1,'first')+di+1;
            end
          
        %end
    end
end


%calculate threshold for each AP
%threshold defined as 5% of the maximum dV/dt
for si = 1:sp_num
    max_dv(si,1) = max(dV_sec(sp_index(si,1):sp_index(si,2)));
end

for si = 1:sp_num
     ap_dv = dV_sec(sp_index(si,1):sp_index(si,2));
     
     th_ind = find(ap_dv>=0.05*nanmean(max_dv),1,'first');
     V_th(si,1) = th_ind+sp_index(si,1)+1;
     V_th(si,2) = V(V_th(si,1));
     
     f_trough(si,1) = th_ind+sp_index(si,2)+1;
     f_trough(si,2) = V(f_trough(si,1));
end

%upstroke/downstroke ratio
for ui = 1:sp_num
    udratio(ui,1) = abs(max(dV_sec(sp_index(ui,1):sp_index(ui,2)))...
    /min(dV_sec(sp_index(ui,1):sp_index(ui,2))));
end
     
%%%%Plotting
if plot_on == 1
    figure(); 
    hold on;
    plot(V,'linewidth',2);
    plot(dV_sec./10,'linewidth',1.5);
    for vi = 1:sp_num
        plot(V_th(vi,1),V_th(vi,2),'o','markersize',8,'markerfacecolor','k');
        plot(f_trough(vi,1),f_trough(vi,2),'o','markersize',8,'markerfacecolor','r');
        hold on
    end
    set(gca,'xlim',[8000 25000]);
end
end