%Cm and Rin data combined from mini and f-I experiments


fp_analyzed_data = 'C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di';

filename = {'chronic_DREADDs_24_5_combined_pp.mat'};

import_data = ...,
    xlsread('C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di\PassiveProperties_mini_fI_combined.xlsx');

Cm_all = import_data(:,1:3);
Rin_all = import_data(:,5:7);