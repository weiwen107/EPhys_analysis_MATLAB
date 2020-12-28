%%%% This script calculates the kinetics and event charge of waveform averages of minis obtained from a
%%%% specific cell. WAVG and wavgs_raw yielded from MINIANALYSIS should be
%%%% loaded into the workspace before running.
%% change these for each run
%location of MINIANALYSIS mat files
fp_mini_analysis = ...,
    'C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di\analyzed_mini_results\';

%rise time cutoff
rise = 'rise_1';

%Name of experiment(s) to be run
experiment = '201029';

%location where kinetic results are saved
fp_kinetics = ...,
    'C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di\waveform_average_kinetics\';

%%%%%plotting settings
figure_on = 0;

%%%%save results
save_results = 1;

%calc_mode = 1: calculate the average waveform across all traces of the
%same cell and then get kinetics of this wave. calc_mode = 2: calculate
%kinetics for each trace and average them across the same cell.
calc_mode = 1; 

%%

cd(strcat(fp_mini_analysis,rise))
load(strcat('MINIANALYSIS_',experiment,'.mat'))

%initiation
riseAve = cell(1,numel(WAVG));
decayAve = cell(1,numel(WAVG));
decayTau = cell(1,numel(WAVG));
decayTot = cell(1,numel(WAVG));
wavg_per_cell = cell(1,numel(WAVG));
wavg_charge = cell(1,numel(WAVG));
coefs = cell(1,numel(WAVG));
gofs = cell(1,numel(WAVG));

%Parameter settings
sample_rate = 5000;
rise_start_pc = 0.1;
rise_end_pc = 0.9;
decay_start_pc = 0.9;
num_pts_for_decayfit = 50;
num_pts_BL = 30;


%%%%Fitting parameters for decay phase
s = fitoptions('Method','NonlinearLeastSquares','Startpoint',[11 -0.4],...
    'Lower',[5 -3],'Upper',[120 -0.04]);
f = fittype('a*exp(b*x)','options',s);


for ei = 1:numel(wavgs_raw)
    for ci = 1:numel(wavgs_raw{1,ei})
        wavg_all_per_cell = [];
        
        if isempty(WAVG{1,ei}{1,ci}) == 1
            continue
        end
        
        if calc_mode == 1
            
            for ti = 1:numel(wavgs_raw{1,ei}{1,ci})
                if isempty(wavgs_raw{1,ei}{1,ci}{1,ti}) == 1
                    continue
                elseif isnan(wavgs_raw{1,ei}{1,ci}{1,ti}) == 1
                    continue
                else
                
                    for mi = 1:size(wavgs_raw{1,ei}{1,ci}{1,ti},1)
                        for ni = 1:size(wavgs_raw{1,ei}{1,ci}{1,ti},2)

                            if wavgs_raw{1,ei}{1,ci}{1,ti}(mi,ni) == 0
                                wavgs_raw{1,ei}{1,ci}{1,ti}(mi,ni) = NaN;
                            end
                        end
                    end

                    wavg_all_per_cell = cat(2,wavgs_raw{1,ei}{1,ci}{1,ti});
                end
            end

            for ii = 1:140%size(wavg_all_per_cell,1)
                wavg_per_cell{1,ei}(ii,ci) = nanmean(wavg_all_per_cell(ii,:));
            end
            
        elseif calc_mode ==2
            wavg_all_per_cell = WAVG{1,ei}{1,ci};
            
            for ii = 1:140%size(wavg_all_per_cell,1)
                wavg_per_cell{1,ei}(ii,ci) = nanmean(wavg_all_per_cell(ii,:));
            end
        end
        
        current_cell_average = wavg_per_cell{1,ei}(:,ci);
        if figure_on == 1
            figure('position',[56 200 1200 490])
            subplot(1,3,1:2)
            plot(0:0.2:0.2*(numel(current_cell_average)-1),current_cell_average)
            title(strcat('Cell_',num2str(ci)),'Interpreter','none')
            ylabel('pA')
            xlabel('ms')
            hold on
        end


        %%%%calculate risetime: use interpolation to get a better resoultion
        
        [amplitude,pk_ind] = min(current_cell_average);
        rise_raw = current_cell_average(pk_ind-14:pk_ind); %use slope to find starting point
        slope = NaN(1,15-2);
        for si = 1:15-2 %calculate slope every 3 points
            slope(si) = (rise_raw(si+2)-rise_raw(si))/2;
        end
        
        interp_start = find(slope < -0.1,1,'first');
        
        interp_end = numel(rise_raw);
        
        interp_rise = interp1(interp_start:interp_end,rise_raw(interp_start:interp_end),interp_start:0.1:interp_end);
        
        rise_start_ind = find(interp_rise <= rise_start_pc*amplitude,1,'first');
        rise_end_ind = find(interp_rise >= rise_end_pc*amplitude,1,'last');
        rise_time = (rise_end_ind-rise_start_ind)/(10*sample_rate);
        interp_time_interval = 1000/(sample_rate/0.1); %in ms
        
        rise_time_stamp = rise_start_ind*interp_time_interval : interp_time_interval : rise_end_ind*interp_time_interval;
        rise_time_stamp = rise_time_stamp + (pk_ind-14-1+interp_start-1)*0.2;
        if figure_on == 1
            plot(rise_time_stamp,interp_rise(rise_start_ind:rise_end_ind),'r','LineWidth',3)
        end
        
        riseAve{1,ei}(ci,1) = rise_time; %in s

        %%%%calculate decaytime
        decay_start_ind = find(current_cell_average(pk_ind:end) >= decay_start_pc*amplitude,1,'first')+pk_ind-1;
        decay_end_ind = decay_start_ind + num_pts_for_decayfit;

        decay_phase = current_cell_average(decay_start_ind:decay_end_ind);
        decayAve{1,ei}(ci,1) = numel(decay_phase)/sample_rate; %in s
        
        if figure_on == 1
            plot(decay_start_ind*0.2:0.2:decay_end_ind*0.2,decay_phase,'b','LineWidth',3)
        end

        %%%%Exponential fitting of decay phase

        [exp_fit,gof] = fit((0:0.2:0.2*num_pts_for_decayfit)',-decay_phase,f);
        coefs{1,ei}(ci,1) = exp_fit.a;
        coefs{1,ei}(ci,2) = exp_fit.b;
        gofs{1,ei}(ci) = gof;

        if figure_on == 1
            subplot(1,3,3)
            plot(0:0.2:0.2*num_pts_for_decayfit,-decay_phase)
            hold on
            plot(exp_fit)
            title('decay phase fitting')
            ylabel('pA')
            xlabel('ms')
            hold off

        end

        decayTau{1,ei}(ci,1) = 1/(-exp_fit.b*1000); %in s
        decayTot{1,ei}(ci,1) = decayTau{1,ei}(ci,1)*log(exp_fit.a); %decay time calculated by exponential fit, in s

        %%%%calculate charge
        wavg_charge{1,ei}(ci,1) = nansum(current_cell_average((pk_ind-14-1+interp_start-1):decay_end_ind))/(-10^12*sample_rate); %in C
        if figure_on == 1
            subplot(1,3,1:2)
            for pi = (pk_ind-14-1+interp_start-1):decay_end_ind
                plot([pi*0.2 pi*0.2],[0 current_cell_average(pi)],'k-')
            end
        end

    end
    filename = strcat(experiment,'_calc_mode',num2str(calc_mode),'_',rise,'.mat');
end

%% save results
if save_results == 1
    cd(strcat(fp_kinetics,rise))
    save(filename,'wavg_per_cell','riseAve','decayAve','decayTau','decayTot','wavg_charge','coefs','gofs')
end
