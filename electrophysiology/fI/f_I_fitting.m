function [mean_MFR,std_MFR,fit_coeffs,fit_results] =  f_I_fitting(current_inj,MFR,plot_on,fitcurve_on)
%%% f_I_fitting(current_inj,MFR,plot_on,fitcurve_on) takes an array of injected
%%% currents and the mean firing rate of each neuron, then calculate the
%%% average firing rate and their standard deviation. It also fits the
%%% linear portion of the data

%%%mean_MFR- the mean MFR at each injected current
%%%std_MFR- the standard deviation of each mean_MFR
%%%fit_coeffs- the 95% confidence intervals of each parameter yielded
%%%fit_results- the goodness of the fit


if nargin < 3
    plot_on = 0;
    fitcurve_on=0;
elseif nargin < 4
    fitcurve_on = 0;
end

mean_MFR = zeros(1,size(MFR,2));
std_MFR = zeros(1,size(MFR,2));


if ~(size(current_inj) == size(MFR,2))
    warning('X and Y should have the same dimension!')
else
    for i = 1:size(MFR,2)
        mean_MFR(i) = mean(MFR(:,i));
        std_MFR(i) = std(MFR(:,i));
    end
    
    
    
    %linear regression after firing has started
    fit_start = find(mean_MFR>0,1,'first');
    
    [f, gof] = fit(current_inj(fit_start:end)',mean_MFR(fit_start:end)','poly1');
    fit_coeffs = confint(f,0.95); % %95 confidence level 
    fit_results = gof;
    interpolated_MFR = f.p1.*current_inj(fit_start:end)+f.p2;
   
    
   if plot_on == 1
    f = figure();
    errorbar(current_inj,mean_MFR,std_MFR,'-o','MarkerSize',10,...
        'MarkerEdgeColor','red', 'MarkerFaceColor','red','LineWidth',2)
    hold on
    
    if fitcurve_on ==1
    plot(current_inj(fit_start:end),interpolated_MFR,'k-','LineWidth',2)
    end
    
    xlabel('Current (pA)');
    xlim([0 current_inj(end)+20]);
    ylabel('Firing Rate (Hz)');
    
    fp = cd;
    total = numel(fp);
    title(strcat('Experiment: ', fp(total-5:end)), 'Interpreter', 'none');
    
    auto_zoom(f);
   end
end
    
