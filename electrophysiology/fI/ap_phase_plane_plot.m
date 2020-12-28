%% import the dV/dt and V of the first AP evoked from each cell and calculate the meta average
save_results = 1;

fp_analyzed_data = 'C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di';
import_file_name = 'firstap_all_cells.xlsx';

save_filename = 'chronic_DREADDs_firstAP_average.mat';

field1 = 'WT_V'; 
field2 = 'WT_CNO_V';
field3 = 'DR_CNO_V';
field4 = 'WT_dV'; 
field5 = 'WT_CNO_dV';
field6 = 'DR_CNO_dV';

data_temp = NaN(150,30);

firstAP_all = struct(field1, data_temp,...
    field2, data_temp,...
    field3, data_temp,...
    field4, data_temp,...
    field5, data_temp,...
    field6, data_temp);

fn = fieldnames(firstAP_all);

for fii = 1:numel(fn)
    firstAP_all.(fn{fii}) =...,
        xlsread(strcat(fp_analyzed_data,'\',import_file_name), fn{fii});
end
        
        
%% calculate average and width at half-height for the first APs

V_th = NaN(30,2);
f_trough = NaN(30,2);

for fi = 1:numel(fn)
    
    curr = firstAP_all.(fn{fi});
    for di = 1:size(curr,1) 
      firstAP_ave.(fn{fi})(di,1) = nanmean(curr(di,(1:size(curr,2))));
    end
    
    if fi<=(numel(fn)/2)
        curr_V = firstAP_all.(fn{fi});
        curr_dV = firstAP_all.(fn{fi+numel(fn)/2});
    
        %width
        for sii = 1:size(curr_V,2)
            [peak, peak_ind] = max(curr_V(:,sii));
            ap_dv = curr_dV(:,sii);
            max_dv = max(ap_dv);

            th_ind = find(ap_dv>=0.05*nanmean(max_dv),1,'first');
            V_th(sii,1) = th_ind;
            V_th(sii,2) = curr_V(V_th(sii,1));

            if isempty(find(ap_dv(peak_ind+10:peak_ind+100) >= 0.01*min(ap_dv(peak_ind+1:peak_ind+100)),1,'first'))
                [ft_val,ft_ind] = min(curr_V(peak_ind+10:peak_ind+80));
                f_trough(sii,1) = ft_ind+peak_ind+1;
            else
                f_trough(sii,1) = find(ap_dv(peak_ind+10:peak_ind+100)>=0.01*min(ap_dv(peak_ind+1:peak_ind+100)),1,'first')+peak_ind+1;
            end

            f_trough(sii,2) = curr_V(f_trough(sii,1),sii);

        %downward halfheight
        halfheight_d = 0.5*(peak - f_trough(sii,2))+f_trough(sii,2);
        halfheight_d_ind = find(curr_V(peak_ind:f_trough(sii,1),sii)<=...
            halfheight_d,1,'first')+peak_ind;

        %extrapolate to the upward side
        data_extrp = curr_V(10:peak_ind,sii);
        [data_extrp, index] = unique(data_extrp);
        %data_extrp = data_extrp';
        ind_extrp = 10:peak_ind;
        ind_extrp = ind_extrp(index)';


        halfheight_u_ind = interp1(data_extrp, ind_extrp, halfheight_d);
        width.(fn{fi})(sii,1) = (halfheight_d_ind - halfheight_u_ind)*0.1;

            %troubleshooting width
    %     figure;
    %     hold on;
    %     plot(curr_V(:,sii),'k')
    %     plot(halfheight_d_ind,halfheight_d,'o','markersize',10,'markerfacecolor','r')
    %     plot(halfheight_u_ind,halfheight_d,'o','markersize',10,'markerfacecolor','b')
        end
    end
end

%% plotting
figure
hold on
for fii = 1:3
    plot(firstAP_ave.(fn{fii}), firstAP_ave.(fn{fii+3}),'LineWidth',2)
    hold on
    %plot(firstAP_ave.WT_dV/10)
end

legend('WT','WT_CNO','DR_CNO','Interpreter','none')
ylabel('dV/dt')
xlabel('mV')

hold off

% figure
% plot(firstAP_ave.WT_V, firstAP_ave.WT_dV)

%% save firstAP average
if save_results == 1
    cd(fp_analyzed_data)
    save(save_filename,'firstAP_all','firstAP_ave','width')
end