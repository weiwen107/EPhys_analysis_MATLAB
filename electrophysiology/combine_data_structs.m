
fields = {'MFR', 'IFR', 'mean_IFR', 'Threshold', 'ADP_ind','Lat',...
    'Udratio','Width_med', 'Width_first','Rheobase', 'Rin', 'Cm','Vm',...,
    'AHP_amp','AHP_dur','AHP_area'}; 

ds1 = adult_WT_CNO;
ds2 = adult_CNO_48h;


for fi = 1:numel(fields)
    if size(ds1.(fields{fi}),2) > 1
        ds_temp.(fields{fi}) = cat(2,ds1.(fields{fi}),ds2.(fields{fi}));
    else
        ds_temp.(fields{fi}) = cat(1,ds1.(fields{fi}),ds2.(fields{fi}));
    end
end

adult_CNO_24_n_48 = ds_temp;
current_inj = (20:20:400)'; %in pA

filename = 'chronic_hm4di_fI_adult_pooled_control_24_48.mat';
save(filename,'adult_CNO_24_n_48','current_inj')