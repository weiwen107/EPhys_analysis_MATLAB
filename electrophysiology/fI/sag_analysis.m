% This script takes hyperpolarization traces obtained under the current clamp
% and calculate the sag amplitude.
% The sag amplitude is defined as the difference between the peak
% hyperpolarized membrane potential and the following steady-state potential 
% during the current injection.
% Steady-state is defined as the region where membrance potential
% is less than 10% of the peak potential.

%% these need to be changed for each new data run

save_file_name = 'sag_200302.mat';
%%%%%% Name of file you would like to have data saved to

%loop through all the folders
experiment = {'200302'};

%filepath where you keep data folders
fp_data = 'C:\Users\schum\Google_Drive\Lab\Data\chronic_hm4di\hyperpolarization\';

%location of the analyzed mini data
fp_analyzed_data = ...,
    'C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di\analyzed_sag';

%save results
save_results = 1;

%plot on
figure_on = 0;

%sample rate
samprate = 10000;

%% data readout and sag amplitude calculation

for jj = 1:1%numel(experiment)
    current_experiment = experiment{jj};
    current_folder = strcat(fp_data,current_experiment);
    cd(current_folder)
    
    all_files = dir;
    file_num = numel(all_files)-2;

    %%% pre-allocate space
    filename = cell(1,file_num);
    cell_id = cell(1,file_num);
    cell_num = NaN(1,file_num);
    aDAT = cell(1,file_num);
    
    peakAmp = cell(1,file_num);
    sagAmp = cell(1,file_num);
    sagAmp_norm = cell(1,file_num);
    sagTau = cell(1,file_num);
    fit_para = cell(1,file_num);
    fit_qual = cell(1,file_num);
    
%     sag_amp_norm_ave = NaN(file_num,1);
%     sag_tau_ave = NaN(file_num,1);
    
    %%% data readout
    for f = 1:file_num
        filename{f} = all_files(f+2).name;
        num_filename = numel(filename{f});

        if ~strfind(filename{f}, '-')
            warning('Not a sag trace!')
        else

            if isnan(str2double(filename{f}(6))) == 1
                cellID = str2double(filename{f}(5));
            else
                cellID = str2double(filename{f}(5:6));
            end
            %cellID = str2double(filename{f}(5));
            cell_num(1,f) = cellID;

            range(1) = str2double(filename{f}(num_filename - 11:num_filename - 8)); %first sweep #
            range(2) = str2double(filename{f}(num_filename - 6:num_filename - 3)); %last sweep #


            %formatting file names
            if range(1) < 10
                start = strcat('000',num2str(range(1)));
            elseif range(1) >= 10 && range(1) < 100
                start = strcat('00',num2str(range(1)));
            elseif range(1) >= 100 && range(1) < 1000
                start = strcat('0',num2str(range(1)));
            else
                start = num2str(range(1));
            end

            if range(2) < 10
                finish = strcat('000',num2str(range(2)));
            elseif range(2) >= 10 && range(2) < 100
                finish = strcat('00',num2str(range(2)));
            elseif range(2) >= 100 && range(2) < 1000
                finish = strcat('0',num2str(range(2)));
            else
                finish = num2str(range(2));
            end


            for tii = range(1):range(2)

                trace_id = tii;
                cell_id{1,cellID}(trace_id,1) = trace_id;

                if tii < 10
                    sname = strcat('000',num2str(tii));
                elseif tii >= 10 && tii < 100
                    sname = strcat('00',num2str(tii));
                elseif tii >= 100 && tii < 1000
                    sname = strcat('0',num2str(tii));
                else
                    sname = num2str(tii);
                end


                %datafile readout
                dtstruct = ws.loadDataFile(filename{f});
                prefix = strcat('dtstruct.sweep_',sname,'.analogScans');
                data = eval(prefix); 
                aDAT{1,cellID}(1:numel(data(:,1)),trace_id) = data(1:numel(data(:,1)),1);

            end
        end
    end

    %%% Replace zeroes in each cell with NaN
    for gi = 1:max(cell_num)
        cell_id{1,gi} = nonzeros(cell_id{1,gi});
        if isempty(cell_id{1,gi}) == 1
            continue
        else
            for ri = 1:cell_id{1,gi}(end,1)
                if ismember(ri,cell_id{1,gi}) == 0
                    aDAT{1,gi}(:,ri) = NaN;
                end
            end
        end
    end
    
    %% sag calculation
    for ci = 1:max(cell_num)
        if isempty(cell_id{1,ci}) == 1
            continue
        else
            for ti = cell_id{1,ci}(1,1):cell_id{1,ci}(end,1)
                
                data = aDAT{1,ci}(8000:22000,ti);
                
%                 figure
%                 plot(data)
%                 hold on
                
                [peak, ind] = min(data);
                
                if ind >= 4000
                    peakAmp{1,ci}(ti,1) = NaN;
                    sagAmp{1,ci}(ti,1) = NaN;
                    sagAmp_norm{1,ci}(ti,1) = NaN;
                    sagTau{1,ci}(ti,1) = NaN;
                    fit_para{1,ci}{ti} = NaN;
                    fit_qual{1,ci}(ti,1) = NaN;
                    
                    disp(strcat('No peak found in Cell',num2str(ci),' Trace',num2str(ti),'!'))
                    
                else
                    ss_est_range = ind+3000 : ind+8000;
                    ss_est_amp = nanmean(data(ss_est_range));
                    %plot(ss_range,data(ss_range))
                    sag_est_amp = ss_est_amp - peak;

                   %%% for the sake of fitting, only the steady-state region right
                   %%% after the rebound should be included (~100 ms)
                    ss_start = find(data(ind:ind+3000) >= peak+sag_est_amp*0.9,1,'first')+ind;
                    ss_end = ss_start + 1000;
                    ss_amp = nanmean(data(ss_start:ss_end));
                    sag_amp = ss_amp - peak;

                    sag_amp_norm = (sag_amp/peak)*100; %normalized sag amplitude (to the peak)

                    % first-order exponential fitting
                    % region should be normalized first
                    fitval = -(data(ind:ss_end)-ss_amp);
                    time = (0:(1/samprate):(length(fitval)/samprate)-(1/samprate)).';

    %                 figure
    %                 plot(fitval)

                    tau_est = (find(fitval < fitval(1)*0.37,1,'first')-1)/samprate;
                    s = fitoptions('Method', 'NonlinearLeastSquares', ...
                        'StartPoint',[fitval(1), tau_est*0.5],...
                        'Lower', [mean(fitval(1:2))*0.9, tau_est*0.01],...
                        'Upper', [fitval(1)*1.1, tau_est*5]);

    %                 s = fitoptions('Method', 'NonlinearLeastSquares', ...
    %                     'StartPoint',[fitval(1), tau_est*0.5, fitval(200), tau_est*1.5],...
    %                     'Lower', [mean(fitval(1:2))*0.9, tau_est*0.01, mean(fitval(1:200))*0.1, tau_est*0.1],...
    %                     'Upper', [fitval(1)*1.1, tau_est*5, fitval(1)*1.25, tau_est*15]);

                    f = fittype('a*exp(-x/b)','options',s);
                    %f = fittype('a*exp(-x/b) + c*exp(-x/d)','options',s);
                    %watchfithappen(time, fitval, f, 10)

                    [exp_fit,gof] = fit(time, fitval, f);

                    cval = coeffvalues(exp_fit);
                    %fit_error = gof.rmse;


                    tau_cal = cval(2); %sag time constant

                    if figure_on == 1
                         figure('position',[56 200 1200 490]); 
                         subplot(1,2,1)
                         plot(data)
                         title(strcat('Cell',num2str(ci),', Trace',num2str(ti)))

                         subplot(1,2,2)
                         plot(exp_fit, time, fitval)
                         title('Fitting')

                    end
                
                
                    if gof.rsquare <=0.8 || gof.rmse >= 0.25
                        peakAmp{1,ci}(ti,1) = NaN;
                        sagAmp{1,ci}(ti,1) = NaN;
                        sagAmp_norm{1,ci}(ti,1) = NaN;
                        sagTau{1,ci}(ti,1) = NaN;
                        fit_para{1,ci}{ti} = gof;
                        fit_qual{1,ci}(ti,1) = 0;
                    else
                        peakAmp{1,ci}(ti,1) = peak;
                        sagAmp{1,ci}(ti,1) = sag_amp;
                        sagAmp_norm{1,ci}(ti,1) = sag_amp_norm;
                        sagTau{1,ci}(ti,1) = tau_cal;
                        fit_para{1,ci}{ti} = gof;
                        fit_qual{1,ci}(ti,1) = 1;
                    end
                end
            end %trace
        end
    end %cell
            
%     sag_amp_norm_ave(ci,1) = nanmean(sagAmp_norm{1,ci}(:,1));
%     sag_tau_ave(ci,1) = nanmean(sagTau{1,ci}(:,1));
            
end %experiment

%% save to file
if save_results == 1
    cd(fp_analyzed_data)
    save(save_file_name, 'peakAmp','sagAmp','sagAmp_norm','sagTau','fit_para','fit_qual')
        %'sag_amp_norm_ave','sag_tau_ave')
end