clear
load JEG_graphics
load ../data/mollusk_db_holo.mat
load ../data/coral_db_holo.mat
addpath ../code/
% analysis parametes
Tmin = 2; Tmax = 7;  % ENSO band boundaries
window = 15;    % window size to compute variance
len_y = 2;   % block bootstrap block length
Nb = 1000;    % block bootstrap # of samples

% order field names
Co = orderfields(C);
Mo = orderfields(M);
fields = fieldnames(Co);

% merge the two databases
Cc = struct2cell(Co);
Mc = struct2cell(Mo);
Sc = cat(3, Cc , Mc);
S  = cell2struct(Sc,fields,1);

% =================================================================================
% query database for suitable records
% =================================================================================
meas   = {S.meas};
O18    = strncmp('\delta^{18}O',meas,12); % find all 18O records
res    = [S.res];
seas   = (res <= 3); % locate seasonally-resolved records
nyears = [S.nyears];

% apply selection if necessary
%Sr = S(O18 & long & seas);
%Sr = S(O18 & seas);
%Sr = S(O18);


%==========================================================
%  ANALYSIS PER SE
%==========================================================
Sr = S; nr = length(Sr);
% placeholder
for r = 1:nr
    Sr(r).seasonal_amp_b = NaN;
    Sr(r).seasonal_amp_ratio_b = NaN;
    Sr(r).enso_var_ratio_b = NaN;
end

% ANNUAL CYCLE AMPLITUDE UNCERTAINTIES
% check that npy is the correct decorrelation length to use. [yes]
ann = find(~isnan([Sr.seasonal_amp]) & isnan([Sr.seasonal_amp_err]));
n_ann = numel(ann); % 49 records 08/07/2014

for a = 1:n_ann
    r = ann(a)
    X = Sr(r).data; t = Sr(r).chron;
    %variance = var(X);
    mRes =  abs(Sr(r).res); % correct for the fact that some records run backwards
    Sr(r).res = mRes; dt = mRes/12; npy = round(1/dt);  mcorr = 1;
    year = floor(t);  month = floor((t - year)*12+mcorr); nm  = length(month);
    tv_reg = [year, month ,repmat(15, [nm 1])];
    tn_reg = datenum(tv_reg);
    period = min(50,Sr(r).nyears);
    
    seas_amp = zeros(Nb,1);
    % resample
    Xb = block_bootstrap(X,npy,Nb);
    for k = 1:Nb
        % remove trend prior to computing seasonal cycle
        [dtd,trend]=spl_dt(t,Xb(:,k),period);
        % compute  seasonal cycle
        [~,Xc] = remove_season2(dtd,tn_reg,1,npy,1); % apply the tropical year option (1)
        seas_amp(k) = range(Xc);
    end
    Sr(r).seasonal_amp_b = seas_amp;
    Sr(r).seasonal_amp_q = quantile(seas_amp,[.025 .25 .5 .75 .975]);
    Sr(r).seasonal_amp_err = std(seas_amp);
end

% rest = setdiff([1:nr],ann);
% %% set others to an empty vector
% for r = rest
%     Sr(r).seasonal_amp_b = [];
% end

% seperate old and modern records
age_mid = [Sr.age_mid]';
So = Sr(age_mid < 1900);   no = length(So);  % old
Sm = Sr(age_mid >= 1900);  nm = length(Sm); % modern
c_lon = [So(:).lon]'; % fix longitudes:
c_lon(c_lon <0) = 360 + c_lon(c_lon <0);
c_lat = [So(:).lat]';
archive = {So(:).archive}';

id_mod = {Sm.modern_id}; % find modern samples
no_ratio = find(~isnan([So.seasonal_amp]) & isnan([So.seasonal_amp_ratio]));

ns = numel(no_ratio);
for s = 1:ns
    r = no_ratio(s);
    fossil_amp_b = So(r).seasonal_amp_b;
    id = strcmp(So(r).modern_id,id_mod);
    modern_amp_b = Sm(id).seasonal_amp_b;
    So(r).seasonal_amp_ratio = So(r).seasonal_amp/Sm(id).seasonal_amp;
    R = fossil_amp_b./modern_amp_b;
    So(r).seasonal_amp_ratio_b = R;
    So(r).seasonal_amp_ratio_q = quantile(R,[.025 .25 .5 .75 .975]);
end


%  VARIANCE RATIO
% ======================
citekey_o = {So.citekey};   % extract bibliographic cite key
citekey_m = {Sm.citekey};   % extract bibliographic cite key
% isolate records that don't need variance calculation
var_calc = find(~strcmp(citekey_o,'Carre_2014'));
r = 1;
while ismember(r,var_calc)
    X = So(r).data; nt = length(X);
    mRes =  So(r).res; % ensure positive resolution   
    len_m = round(12/mRes*len_y);
    % apply wavelet filter
    [Xe,vra] = wavelet_filter(X,mRes/12,Tmin,Tmax);
    % extract enso band variance
    So(r).enso_var = var(Xe);
    
    % estimate sampling distribution via bootstrap
    Xbe = block_bootstrap(Xe,len_m,Nb); 
    
    %  modern analog
    id = strcmp(So(r).modern_id,id_mod);
    Xm = Sm(id).data;
    [Xme,vra] = wavelet_filter(Xm,mRes/12,Tmin,Tmax);
    Xmbe = block_bootstrap(Xme,len_m,Nb);
    
    Vb = var(Xbe); Vmb = var(Xmbe);
    So(r).enso_var_q = quantile(Vb,[.025 .25 .5 .75 .975]);
    So(r).enso_var_err = std(Vb);
    
    %  variance ratio
    R = Vb./Vmb;
    So(r).enso_var_ratio =  var(Xe)/var(Xme);
    So(r).enso_var_ratio_q = quantile(R,[.025 .25 .5 .75 .975]); 
    So(r).enso_var_ratio_err = std(R); 
    So(r).enso_var_ratio_b = R;
    r = r+1
end

% compute ratios for modern samples
noCarre = find(~strcmp(citekey_m,'Carre_2013'));
for r = noCarre % skip the Carré et al samples
    X = Sm(r).data; nt = length(X);
    t = Sm(r).chron;  
    mRes =  abs(Sm(r).res); % ensure positive resolution
    Sm(r).res = mRes;
    
    Sm(r).seasonal_amp_ratio = 1;
    R = Sm(r).seasonal_amp_b/Sm(r).seasonal_amp;
    Sm(r).seasonal_amp_ratio_q =  quantile(R,[.025 .25 .5 .75 .975]);
    Sm(r).seasonal_amp_err = std(R);
    
    %%  variance ratio
    len_m = round(12/mRes*len_y);
    % apply wavelet filter
    [Xe,vra] = wavelet_filter(X,mRes/12,Tmin,Tmax);
    % extract enso band variance
    Sm(r).enso_var = var(Xe);
    % estimate sampling distribution via bootstrap
    Xbe = block_bootstrap(Xe,len_m,Nb);
    Vb = var(Xbe);
    Sm(r).enso_var_q = quantile(Vb,[.025 .25 .5 .75 .975]);
    Sm(r).enso_var_err = std(Vb);
    
    R = Vb/var(Xe);
    Sm(r).enso_var_ratio =  1;
    Sm(r).enso_var_ratio_b = R;
    Sm(r).enso_var_ratio_q =  quantile(R,[.025 .25 .5 .75 .975]);;
    Sm(r).enso_var_ratio_err = std(R);
end

% merge old and new into Sr
Mc = struct2cell(Sm);
Oc = struct2cell(So);
Sc = cat(3, Oc , Mc);
fields = fieldnames(Sm);
Sr  = cell2struct(Sc,fields,1);
citekey = {Sr.citekey};
% Swap records around to put matthieu's records last
Carre = find(strcmp(citekey,'Carre_2013')|strcmp(citekey,'Carre_2014'));
%43    44    45    46    47    48    62
idx = [1:nr];
idx(43:48) = 54:59; idx(54:59) = 43:48;
Sw = Sr(idx); cite_sw = {Sw.citekey};
Sr = Sw;

% export and save
save '../data/proxy_db_holo_synthesis.mat' So Sm Sr

