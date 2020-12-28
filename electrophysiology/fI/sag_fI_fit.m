%% save file
save_file_name = 'fI_sag_data_adult.mat';
%%%%%% Name of file you would like to have data saved to

% where to save the analyzed data
fp_analyzed_data = ...,
    'C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di\fI_data_by_groups';

%save results
save_results = 0;

%% conditions
condn = {'WT_CNO','DR_CNO'};

%vari = {'MFR','IFR','sag_amp','sag_tau'};

%% data arrangement
%datatemp = cell(1,numel(vari));

for gi = 1:numel(condn)
    
    datatemp1 = eval(strcat(condn{gi},'.sag_amp'));
    datatemp2 = eval(strcat(condn{gi},'.sag_tau'));
    datatemp3 = eval(strcat(condn{1,gi},'_area.MFR'));
    datatemp4 = eval(strcat(condn{1,gi},'_area.IFR'));
    
    [~,cell_indx] =count_non_nan(datatemp1);
    counter = 0;
    Atemp = NaN(size(datatemp1,1),6);
    
    for ci = 1:size(datatemp1,1)
        if ismember(ci,cell_indx)
            counter = counter+1;
            Atemp(counter,1) = datatemp3(ci,1); %col 1: MFR area
            Atemp(counter,2) = datatemp4(ci,1); %col 2: IFR area
            Atemp(counter,3) = datatemp1(ci,1); %col 3: sag_amp
            Atemp(counter,4) = datatemp2(ci,1); %col 4:sag_tau
            Atemp(counter,5) = datatemp1(ci,1)*datatemp2(ci,1)*1000; %area: mv*ms
            Atemp(1,6) = counter; %number of non NaN elements (for fitting purpose)
        else
            continue
        end
    end

    sagdata.(condn{gi}) = Atemp;
    
end
 %% plotting

 scale = 1;
 X = cat(1,sagdata.WT_CNO(1:sagdata.WT_CNO(1,6),3),sagdata.DR_CNO(1:sagdata.DR_CNO(1,6),3)).*-1;
 Y = cat(1,sagdata.WT_CNO(1:sagdata.WT_CNO(1,6),4),sagdata.DR_CNO(1:sagdata.DR_CNO(1,6),4)).*scale;
 
%  [P,S,mu] = polyfit(X,Y,1);
%  Xfit = linspace(min(X),max(X),100);
%  Yfit = polyval(P,Xfit);
[f,gof] = fit(X,Y,'poly1'); 
p = corr(X,Y);
 figure()
 plot(X,Y,'o')
 hold on
 plot(f)
 title('sag amp vs. sag tau')
 
 %% save files
 if save_results == 1
     cd(fp_analyzed_data)
     save(save_file_name,'sagdata')
 end