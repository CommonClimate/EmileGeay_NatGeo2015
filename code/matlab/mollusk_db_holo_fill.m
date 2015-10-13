% This script fills in the data structure M for mollusq 
clear all; close all; 
fields={'archive','meas','sample_id','genus','species','chron','data','anom','year_i','year_f','res','lat','lon','site','reference','citekey','age_mid','age_std'};
Nf=length(fields);  
% TRIDACNA
% from Tridacna_data_set.xls (sent by M. Elliott 19/03/2013)
% ==================                                     Age   Uncertainty
%HU-04-MT7	Modern	Tridacna gigas	Inner			               0			0
%HU-04-T58	Holocene	Tridacna gigas	Hinge	6		    ?     8.09	0.15	   0.1		0.085		
%HU-04-T66	Holocene	Tridacna gigas	Inner	20	      6480	6.54	0.38	   0.6		0
%HU-04-T73	Holocene	Tridacna gigas	Hinge	Ex situ	7129	7.28	0.12		0.1		0.043		
%HU-04-T75	Holocene	Tridacna gigas	Hinge	Ex situ	8592	8.76	0.46		0.2		0.128		
% further precisions from Robin Driscoll (27/03/2013):
% 1) The location is 147.5E and 6.5S
%2) The modern sample is from 1986-2002
%3) The uncertainty refers to 2-sigma error
%  Interpetation: d18O signal is about 50% SST, 50% isotopic dilution due
%  to rainfall. Both effects are synergistic for annual cycle and ENSO.
% Note: no way to get amplitude of annual cycle


% ==============================================
%   MT7 sample (modern)
% ==============================================
n = 1;
M(n).archive   = 'mollusk';  M(n).meas    = '\delta^{18}O';
M(n).genus     = 'Tridacna'; M(n).species = 'gigas'; 
M(n).reference = 'Welsh et al [2011]'; 
M(n).citekey   = 'Welsh_2011';
M(n).site   = 'Huon Peninsula, PNG';
M(n).lon  = 147.5; M(n).lat  = -6.5;
%  fill in specific information
M(n).sample_id = 'HU04MT7';
M(n).modern_id = 'HU04MT7';
[MT7,TXT,RAW] = xlsread('../data/HU04MT7.xls');  % about 16 years, no trend
M(n).year_i = 1986; M(n).age_std = 0; % no uncertainty within machine precision. 
M(n).chron = M(n).year_i + MT7(:,1); t_raw = M(n).chron; 
[X,t,Xr,tr,mr,yr,npy,mRes,dt] = coral_interp(MT7(:,1),MT7(:,2));
[Xa,Xc] = coral_remove_seas(t,X,tr,Xr,mr,yr,npy,M(n).site);

M(n).nyears = ceil(max(M(n).chron)) - floor(min(M(n).chron))+1;
M(n).year_f = ceil(max(M(n).chron));
M(n).age_mid = median(M(n).chron);
M(n).anom = Xa;
M(n).data = M(n).anom;
M(n).res = 12*mode(diff(M(n).chron));
M(n).seasonal_amp = range(Xc);

% ==============================================
%   T58 sample
% ==============================================
n = n+1;
M(n) = M(n-1); % fill in generalities
% fill in specifics
M(n).reference = 'Driscoll [2013]'; 
M(n).citekey   = 'Driscoll_2013';

% REDO with detrended data 

M(n).sample_id = 'HU04T58'; 
M(n).modern_id = 'HU04MT7';
T58 = csvread('../data/HU04T58_new.csv');  % ~22 years, huge trend
[X,t,Xr,tr,mr,yr,npy,mRes,dt] = coral_interp(T58(:,1),T58(:,2));
[Xa,Xc] = coral_remove_seas(t,X,tr,Xr,mr,yr,npy,M(n).site);

M(n).year_i = 1950 - 8.09*1000; 
M(n).age_std = 150/2;  % divide by two to obtain 1-sigma error
M(n).chron = M(n).year_i + t; 
M(n).nyears = ceil(max(M(n).chron)) - floor(min(M(n).chron))+1;
M(n).year_f = ceil(max(M(n).chron));
M(n).age_mid = median(M(n).chron);
M(n).anom = Xa;
M(n).data = M(n).anom;
M(n).res = 12*mode(diff(M(n).chron));
M(n).seasonal_amp = range(Xc);

% ==============================================
%   T66 sample : TOO SHORT
% ==============================================
% n = n+1;
% M(n) = M(n-1); % fill in generalities
% % fill in specifics
% M(n).sample_id = 'HU04T66';
% M(n).modern_id = 'HU04MT7';
% T66           = csvread('../data/HU04T66.csv',1,0,'A2..B79');% ~7 years, no trend
% M(n).year_i = 1950 - 6.54*1000;  
% M(n).age_std = 380/2;  % divide by two to obtain 1-sigma error
% M(n).chron = M(n).year_i + T66(:,1); 
% M(n).nyears = ceil(max(M(n).chron)) - floor(min(M(n).chron))+1;
% M(n).year_f = ceil(max(M(n).chron));
% M(n).age_mid = median(M(n).chron);
% M(n).anom = T66(:,2);
% M(n).data = M(n).anom;
% M(n).res = 12*mode(diff(M(n).chron));
% M(n).seasonal_amp = NaN;

% ==============================================
%   T75 sample
% ==============================================
n = n+1;
M(n) = M(n-1); % fill in generalities
% fill in specifics
M(n).sample_id = 'HU04T75';
M(n).modern_id = 'HU04MT7';
T75 = csvread('../data/HU04T75_new.csv');  % ~22 years, huge trend
[X,t,Xr,tr,mr,yr,npy,mRes,dt] = coral_interp(T75(:,1),T75(:,2));
[Xa,Xc] = coral_remove_seas(t,X,tr,Xr,mr,yr,npy,M(n).site);

M(n).year_i = 1950 - 8.76*1000; 
M(n).age_std = 460/2;  % divide by two to obtain 1-sigma error
M(n).chron = M(n).year_i + t; 
M(n).nyears = ceil(max(M(n).chron)) - floor(min(M(n).chron))+1;
M(n).year_f = ceil(max(M(n).chron));
M(n).age_mid = median(M(n).chron);
M(n).anom = Xa;
M(n).data = M(n).anom;
M(n).res = 12*mode(diff(M(n).chron));
M(n).seasonal_amp = range(Xc);

% [T75,TXT,RAW] = xlsread('../data/HU04T75.xls');  % ~13 years, variance fanning out
% M(n).year_i = 1950 - 8.76*1000;  
% M(n).age_std = 460/2;  % divide by two to obtain 1-sigma error
% M(n).chron = M(n).year_i + T75(:,1); 
% M(n).nyears = ceil(max(M(n).chron)) - floor(min(M(n).chron))+1;
% M(n).year_f = ceil(max(M(n).chron));
% M(n).anom = T75(:,2);
% M(n).age_mid = median(M(n).chron);
% M(n).data = M(n).anom;
% M(n).res = 12*mode(diff(M(n).chron));
% M(n).seasonal_amp = NaN;
% t_raw = M(n).chron;
% X_raw = M(n).data;
% [X,t,Xr,tr,mr,yr,npy,mRes,dt] = coral_interp(t_raw,X_raw);
% % remove seasonal cycle, if any
% [Xa,Xc] = coral_remove_seas(t,X,tr,Xr,mr,yr,npy,M(n).site);

% ==============================================
%   T73 sample
% ==============================================
n = n+1;
M(n) = M(n-1); % fill in generalities
% fill in specifics
M(n).sample_id = 'HU04T73';
M(n).modern_id = 'HU04MT7';
T73 = csvread('../data/HU04T73_new.csv');  % ~22 years, huge trend
[X,t,Xr,tr,mr,yr,npy,mRes,dt] = coral_interp(T73(:,1),T73(:,2));
[Xa,Xc] = coral_remove_seas(t,X,tr,Xr,mr,yr,npy,M(n).site);

M(n).year_i = 1950 - 7.28*1000; 
M(n).age_std = 120/2;  % divide by two to obtain 1-sigma error 
M(n).chron = M(n).year_i + t; 
M(n).nyears = ceil(max(M(n).chron)) - floor(min(M(n).chron))+1;
M(n).year_f = ceil(max(M(n).chron));
M(n).age_mid = median(M(n).chron);
M(n).anom = Xa;
M(n).data = M(n).anom;
M(n).res = 12*mode(diff(M(n).chron));
M(n).seasonal_amp = range(Xc);

% [T73,TXT,RAW] = xlsread('../data/HU04T73.xls');  
% M(n).year_i = 1950 - 7.28*1000;   	
% M(n).age_std = 120/2;  % divide by two to obtain 1-sigma error
% M(n).chron = M(n).year_i + T73(:,1); 
% M(n).nyears = ceil(max(M(n).chron)) - floor(min(M(n).chron))+1;
% M(n).year_f = ceil(max(M(n).chron));
% M(n).anom = T73(:,2);
% M(n).age_mid = median(M(n).chron);
% M(n).data = M(n).anom;
% M(n).res = 12*mode(diff(M(n).chron));
% M(n).seasonal_amp = NaN;

% ======================================
%  apply generic rules to all records
% =====================================
% define analysis parameters
%Tmin = 2; Tmax = 7;  len = 12*2;
%id_mod = {M.modern_id}; mRes = 1;
%Nb = 200;
% placeholders
for j = 1:n
    % enso var 
    M(j).enso_var     = NaN;
    M(j).enso_var_q   = NaN;
    M(j).enso_var_err = NaN;
    % seasonal amplitude
    M(j).seasonal_amp_q   = NaN;
    M(j).seasonal_amp_err = NaN;
    % enso var ratio
    M(j).enso_var_ratio     = NaN;
    M(j).enso_var_ratio_q   = NaN;
    M(j).enso_var_ratio_err = NaN;
    % seasonal amplitude ratio
    M(j).seasonal_amp_ratio     = NaN;
    M(j).seasonal_amp_ratio_q   = NaN;
    M(j).seasonal_amp_ratio_err = NaN;
end


% =================================================================
%  MESODESMA DONACIUM from the Peruvian coast
%  sent by Matthieu Carré on 01/04/13 (no joke!)
%  file Peru_ENSO-Var_dataset.xls
% ==================================================================
fname = '/Users/jeg/Dropbox/PalaeoVar_Working_Folder/MH_synthesis/data/Peru_ENSO-Var_dataset.xls';
[Peru,TXT,RAW] = xlsread(fname,'reduced-grouped'); 
%[Peru-r,TXT,RAW] = xlsread('../data/Peru_ENSO-Var_dataset.xls','reduced-grouped'); 
sigma = [Peru(:,4) Peru(:,6)];
load('../data/Quantiles_Carre_orig.mat');
load('../data/Quantiles_Carre_ratios.mat');

% Modern sample
% ==============
n = n+1; np = n;
M(n).archive   = 'mollusk';  M(n).meas    = '\delta^{18}O';
M(n).genus     = 'Mesodesma'; M(n).species = 'donacium'; 
M(n).reference = 'Carré et al., 2013'; 
M(n).citekey   = 'Carre_2013';
M(n).site   = 'Peruvian Coast';
M(n).lon  = Peru(1,2); 
M(n).lat  = Peru(1,1);        
M(n).sample_id = 'ICA-mod + Llostay';
M(n).modern_id = 'ICA-mod + Llostay';
M(n).nyears       = Peru(1,15);  
M(n).year_i = 1950 - (Peru(1,3)+M(n).nyears/2);
M(n).year_f = 1950 - (Peru(1,5)-M(n).nyears/2);
M(n).age_mid =  1950 - Peru(1,7);
M(n).age_std = max(sigma(1,:))/2; % modern data is certain. 
M(n).res    = 1; % monthly resolution
M(n).enso_var     = Qvar; % 
M(n).enso_var_err = Peru(1,21); %
M(n).enso_var_ratio     = Peru(1,22);
M(n).enso_var_ratio_err = Peru(1,23);
M(n).seasonal_amp              = 3.44*0.23;
M(n).seasonal_amp_err          = 0.16*0.23;
M(n).seasonal_amp_ratio        = 1;
M(n).seasonal_amp_ratio_err    = 0;

% Chala+IS2+Lomas
% ==============
n = n+1;
M(n).archive   = 'mollusk';  M(n).meas    = '\delta^{18}O';
M(n).genus     = 'Mesodesma'; M(n).species = 'donacium'; 
M(n).reference = 'Carré et al., 2014'; 
M(n).citekey   = 'Carre_2014';
M(n).site   = 'Peruvian Coast';
M(n).lon  = Peru(2,2); 
M(n).lat  = Peru(2,1);        
M(n).sample_id = 'Chala+IS2+Lomas';
M(n).modern_id = 'ICA-mod + Llostay';
M(n).nyears       = Peru(2,15);
M(n).year_i = 1950 - (Peru(2,3)+M(n).nyears/2);
M(n).year_f = 1950 - (Peru(2,5)-M(n).nyears/2);
M(n).age_mid =  1950 - Peru(2,7);
M(n).age_std = sqrt(harmmean(sigma(2,:).^2));
M(n).res    = 1; % monthly resolution
M(n).enso_var     = Peru(2,18);
M(n).enso_var_err = Peru(2,21);
M(n).enso_var_ratio     = Peru(2,22);
M(n).enso_var_ratio_q   = Qvar2(2,:);
M(n).enso_var_ratio_err = Peru(2,23);
M(n).seasonal_amp              = 2.70*0.23;
M(n).seasonal_amp_err          = 0.11*0.23;
M(n).seasonal_amp_ratio        = 0.78;
M(n).seasonal_amp_ratio_err   = 0.05;

% Asia 3+5+Ancon
% ==============
n = n+1;
M(n).archive   = 'mollusk';  M(n).meas    = '\delta^{18}O';
M(n).genus     = 'Mesodesma'; M(n).species = 'donacium'; 
M(n).reference = 'Carré et al., 2014'; 
M(n).citekey   = 'Carre_2014';
M(n).site   = 'Peruvian Coast'; 
M(n).lon  = Peru(3,2); 
M(n).lat  = Peru(3,1);        
M(n).sample_id = 'Asia 3+5+Ancon';
M(n).modern_id = 'ICA-mod + Llostay';
M(n).nyears       = Peru(3,15);
M(n).year_i = 1950 - (Peru(3,3)+M(n).nyears/2);
M(n).year_f = 1950 - (Peru(3,5)-M(n).nyears/2);
M(n).age_mid =  1950 - Peru(3,7);
M(n).age_std = sqrt(harmmean(sigma(3,:).^2));
M(n).res    = 1; % monthly resolution
M(n).enso_var     = Peru(3,18);
M(n).enso_var_err = Peru(3,21);
M(n).enso_var_ratio     = Peru(3,22);
M(n).enso_var_ratio_err = Peru(3,23);
M(n).seasonal_amp              = 3.43*0.23;
M(n).seasonal_amp_err          = 0.16*0.23;
M(n).seasonal_amp_ratio        = 0.99;
M(n).seasonal_amp_ratio_err   = 0.07;

% QLB-H4
% ==============
n = n+1;
M(n).archive   = 'mollusk';  M(n).meas    = '\delta^{18}O';
M(n).genus     = 'Mesodesma'; M(n).species = 'donacium'; 
M(n).reference = 'Carré et al., 2014'; 
M(n).citekey   = 'Carre_2014';
M(n).site   = 'Peruvian Coast'; 
M(n).lon  = Peru(4,2); 
M(n).lat  = Peru(4,1);        
M(n).sample_id = 'QLB-H4';
M(n).modern_id = 'ICA-mod + Llostay';
M(n).nyears       = Peru(4,15);
M(n).year_i = 1950 - (Peru(4,3)+M(n).nyears/2);
M(n).year_f = 1950 - (Peru(4,5)-M(n).nyears/2);
M(n).age_mid =  1950 - Peru(4,7);
M(n).age_std = sqrt(harmmean(sigma(4,:).^2));
M(n).res    = 1; % monthly resolution
M(n).enso_var     = Peru(4,18);
M(n).enso_var_err = Peru(4,21);
M(n).enso_var_ratio     = Peru(4,22);
M(n).enso_var_ratio_err = Peru(4,23);
M(n).seasonal_amp              = 2.25*0.23;
M(n).seasonal_amp_err          = 0.24*0.23;
M(n).seasonal_amp_ratio        = 0.65;
M(n).seasonal_amp_ratio_err   = 0.08;

% ICA-IN+QLBN2+3
% ==============
n = n+1;
M(n).archive   = 'mollusk';  M(n).meas    = '\delta^{18}O';
M(n).genus     = 'Mesodesma'; M(n).species = 'donacium'; 
M(n).reference = 'Carré et al., 2014'; 
M(n).citekey   = 'Carre_2014';
M(n).site   = 'Peruvian Coast'; 
M(n).lon  = Peru(5,2); 
M(n).lat  = Peru(5,1);        
M(n).sample_id = 'ICA-IN+QLBN2+3';
M(n).modern_id = 'ICA-mod + Llostay';
M(n).nyears       = Peru(5,15);
M(n).year_i = 1950 - (Peru(5,3)+M(n).nyears/2);
M(n).year_f = 1950 - (Peru(5,5)-M(n).nyears/2);
M(n).age_mid =  1950 - Peru(5,7);
M(n).age_std = sqrt(harmmean(sigma(5,:).^2));
M(n).res    = 1; % monthly resolution
M(n).enso_var     = Peru(5,18);
M(n).enso_var_err = Peru(5,21);
M(n).enso_var_ratio     = Peru(5,22);
M(n).enso_var_ratio_err = Peru(5,23);
M(n).seasonal_amp              = 2.88*0.23;
M(n).seasonal_amp_err          = 0.16*0.23;
M(n).seasonal_amp_ratio        = 0.84;
M(n).seasonal_amp_ratio_err   = 0.06;

% QLB-N4 + QLBN5
% ==============
n = n+1;
M(n).archive   = 'mollusk';  M(n).meas    = '\delta^{18}O';
M(n).genus     = 'Mesodesma'; M(n).species = 'donacium'; 
M(n).reference = 'Carré et al., 2014'; 
M(n).citekey   = 'Carre_2014';
M(n).site   = 'Peruvian Coast'; 
M(n).lon  = Peru(6,2); 
M(n).lat  = Peru(6,1);        
M(n).sample_id = 'QLB-N4 + QLBN5';
M(n).modern_id = 'ICA-mod + Llostay';
M(n).nyears       = Peru(6,15);
M(n).year_i = 1950 - (Peru(6,3)+M(n).nyears/2);
M(n).year_f = 1950 - (Peru(6,5)-M(n).nyears/2);
M(n).age_mid =  1950 - Peru(6,7);
M(n).age_std = sqrt(harmmean(sigma(6,:).^2));
M(n).res    = 1; % monthly resolution
M(n).enso_var     = Peru(6,18);
M(n).enso_var_err = Peru(6,21);
M(n).enso_var_ratio     = Peru(6,22);
M(n).enso_var_ratio_err = Peru(6,23);
M(n).seasonal_amp              = 2.80*0.23;
M(n).seasonal_amp_err          = 0.20*0.23;
M(n).seasonal_amp_ratio        = 0.81;
M(n).seasonal_amp_ratio_err   = 0.07;

% QLB-N6 + QLBN7
% ==============
n = n+1;
M(n).archive   = 'mollusk';  M(n).meas    = '\delta^{18}O';
M(n).genus     = 'Mesodesma'; M(n).species = 'donacium'; 
M(n).reference = 'Carré et al., 2014'; 
M(n).citekey   = 'Carre_2014';
M(n).site   = 'Peruvian Coast'; 
M(n).lon  = Peru(7,2); 
M(n).lat  = Peru(7,1);        
M(n).sample_id = 'QLB-N6 + QLBN7';
M(n).modern_id = 'ICA-mod + Llostay';
M(n).nyears       = Peru(7,15);
M(n).year_i = 1950 - (Peru(7,3)+M(n).nyears/2);
M(n).year_f = 1950 - (Peru(7,5)-M(n).nyears/2);
M(n).age_mid =  1950 - Peru(7,7);
M(n).age_std = sqrt(harmmean(sigma(7,:).^2));
M(n).res    = 1; % monthly resolution
M(n).enso_var     = Peru(7,18);
M(n).enso_var_err = Peru(7,21);
M(n).enso_var_ratio     = Peru(7,22);
M(n).enso_var_ratio_err = Peru(7,23);
M(n).seasonal_amp              = 2.43*0.23;
M(n).seasonal_amp_err          = 0.14*0.23;
M(n).seasonal_amp_ratio        = 0.71;
M(n).seasonal_amp_ratio_err    = 0.05;

% assign quantiles using Matthieu's data
for j = np:n
    jp = j - np +1;
    M(j).enso_var_q = Qvar(jp,:);
    M(j).enso_var_ratio_q = Qvar2(jp,:);
    M(j).seasonal_amp_q = Qseason(jp,:);
    M(j).seasonal_amp_ratio_q = Qseason2(jp,:);
end

% save it all
save '../data/mollusk_db_holo.mat' M



