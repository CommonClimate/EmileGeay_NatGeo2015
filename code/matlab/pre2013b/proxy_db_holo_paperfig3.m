clear, load JEG_graphics
load  '../../../data/obs/proxy_db_holo_synthesis.mat' 

%  Select seasonally resolved records
res = [Sr.res];
Ss  = Sr(res <=3); 
ns = length(Ss);
archive = {Ss.archive};
site    = {Ss.site};
ref     = {Ss.reference};
cite    = {Ss.citekey};
c_lon = [Ss.lon]'; % fix longitudes:
c_lon(c_lon <0) = 360 + c_lon(c_lon <0);
c_lat = [Ss.lat]';

% define regions
reg_ind{1} = find(c_lon <=180); 
reg_ind{2} = find(c_lon > 180 & c_lon < 240);
reg_ind{3} = find(c_lon > 240);

% =============================================
% Link between ENSO-band variance & seasonality
% ==============================================
%err = [So.enso_var_ratio_err]; nyears = [So.nyears];
%res  = [So.res]; nr = length(res);
% Define shortlist
%shortlist = find(err < 1 & res <= 3 & nyears>=20);  % 37 records
%shortlist  = [1:nr]; 
%Ss = So(shortlist); ns = length(Ss);

% assign plotting parameters
for r = 1:ns
   if strfind(site{r},'Vanuatu')
      Ss(r).col = rgb('SaddleBrown'); Vanu = r;
   elseif sum(strncmp(site{r},{'Baye','Surp'},4))>0
      Ss(r).col = rgb('Crimson');     Thierry = r;
   elseif strncmp(ref{r},{'Cobb'},4)
      Ss(r).col = rgb('LimeGreen');   Kim = r;
   elseif strncmp(ref{r},{'Wood'},4)
      Ss(r).col = rgb('LawnGreen');   Woodroffe = r;
   elseif strfind(ref{r},'Gregor')
      Ss(r).col  = rgb('DarkViolet');     Helen = r;
   elseif strfind(ref{r},'Tudhope')
      Ss(r).col  = rgb('DarkOrange'); Sandy = r;
   elseif strfind(ref{r},'Carr')
      Ss(r).col  = rgb('SteelBlue');   Matthieu = r;
   elseif strncmp(ref{r},'Welsh',5) | strncmp(ref{r},'Driscoll',8)
      Ss(r).col  = rgb('DeepSkyBlue'); Robin = r;  %rgb('DeepSkyBlue');
   end
    % assign plotting parameters
   if strcmp(archive{r},'coral')
      mark(r)='o';   sz(r) = 12;
   elseif strcmp(archive{r},'mollusk');
      mark(r)='p';   sz(r) = 16;
   end
end
%
sz2 = round(0.8*sz);
%
col = {Ss.col}; 
lbl{1} = 'Duprey et al [2012]';
lbl{2} = 'this study';
lbl{3} = 'Cobb et al. [2013]';
lbl{4} = 'Woodroffe et al. [2003]';
lbl{5} = 'McGregor et al. [2004,2013]';
lbl{6} = 'Tudhope et al. [2001]';
lbl{7} = 'Carré et al. [2014]';
lbl{8} = 'Driscoll et al. [2014]';
%
% regression variables
Xo = [Ss.seasonal_amp_ratio];  Yo = [Ss.enso_var_ratio];
Nb = 1000; % number of bootstrap samples

%  load data from proxies
% ================================================
for s = 1:ns
    width  = Ss(s).seasonal_amp_ratio_q(4) - Ss(s).seasonal_amp_ratio_q(2);
    seasonal_amp_err(s) = width;
    height = Ss(s).enso_var_ratio_q(4) - Ss(s).enso_var_ratio_q(2);
    enso_var_err(s)     = height;
end

% rescale variables by their expected errors
Xos = Xo./seasonal_amp_err;
Yos = Yo./enso_var_err;
% fit TLS model
[Po,Eo,Xo_hats,Yo_hats] = fitTLS(Xos, Yos);
% boostrap estimate of regression coefficients
Po_boot = bootstrp(Nb, @fitTLS,Xos,Yos);
Po_se = std(Po_boot(:,1)); % slope comes first, then offset. 

% generate family of lines
del = 0.5;
xo = [0:del:20]; nxo = numel(xo);
Yo_hats = polyval(Po,xo);
Lo = zeros(Nb,nxo); %matrix of regression lines
for b = 1:Nb
    Lo(b,:) = polyval(Po_boot(b,:),xo);
end
yo = [-2:del:8]; nyo = numel(yo);
Lo_dens = zeros(nxo,nyo);
for j = 1:nxo
    Lo_dens(j,:)= ksdensity(Lo(:,j),yo);
end
Loq = quantile(Lo,[0.025 0.975],1);

%colormap(brewermap(21,'Purples'))
%  load data from proxies PMIP3 GCMs
% ==========================================
load ../../../data/pmip3/PMIP3_ratios_distrib3.mat
nm = length(models);
for j = 1:3
    Xm(1:nm,j)          = seas_amp_ratio.quant.PIHT(:,j,3);
    %Xm_wdn(1:nm,j)      = seas_amp_ratio.quant.PIHT(:,j,1);
    %Xm_w(1:nm,j)        = seas_amp_ratio.quant.PIHT(:,j,5)-seas_amp_ratio.quant.PIHT(:,j,1);
    Xm_iqr(1:nm,j)      = seas_amp_ratio.quant.PIHT(:,j,4)-seas_amp_ratio.quant.PIHT(:,j,2);
    Xm(nm+1:2*nm,j)     = seas_amp_ratio.quant.MHHT(:,j,3);
    %Xm_wdn(nm+1:2*nm,j) = seas_amp_ratio.quant.MHHT(:,j,1);
    %Xm_w(nm+1:2*nm,j)   = seas_amp_ratio.quant.MHHT(:,j,5)-seas_amp_ratio.quant.MHHT(:,j,1);
    Xm_iqr(nm+1:2*nm,j) = seas_amp_ratio.quant.MHHT(:,j,4)-seas_amp_ratio.quant.MHHT(:,j,2);
    %
    Ym(1:nm,j)          = enso_var_ratio.quant.PIHT(:,j,3);
    %Ym_wdn(1:nm,j)      = enso_var_ratio.quant.PIHT(:,j,2);
    %Ym_w(1:nm,j)        = enso_var_ratio.quant.PIHT(:,j,5)-enso_var_ratio.quant.PIHT(:,j,1);
    Ym_iqr(1:nm,j)      = enso_var_ratio.quant.PIHT(:,j,4)-enso_var_ratio.quant.PIHT(:,j,2);
    Ym(nm+1:2*nm,j)     = enso_var_ratio.quant.MHHT(:,j,3);
    %Ym_wdn(nm+1:2*nm,j) = enso_var_ratio.quant.MHHT(:,j,2);
    %Ym_w(nm+1:2*nm,j)   = enso_var_ratio.quant.MHHT(:,j,5)-enso_var_ratio.quant.MHHT(:,j,1);
    Ym_iqr(nm+1:2*nm,j) = enso_var_ratio.quant.MHHT(:,j,4)-enso_var_ratio.quant.MHHT(:,j,2);
end

%  graphical definitions
colm = flipud(hsv(nm)); colm = repmat(colm,[2 1]); 

% PLOT RAW OBS
win = [0 3 0 2];
fig('ENSO_var vs seasonality, RAW'), clf
subplot(2,1,1)
hold on
for s = 1:ns
    w  = seasonal_amp_err(s);
    he = enso_var_err(s);
    hr(s)  = rectangle('Curvature', [1 1],'Position',[Xo(s)-w/2 Yo(s)-he/2 w he],'Edgecolor',col{s},'Facecolor','none','LineWidth',1);    
    hp(s)  = plot(Xo(s),Yo(s),'marker',mark(s),'MarkerFaceColor',col{s},'MarkerEdgeColor','k','MarkerSize',sz2(s));
end
fancyplot_deco('Proxy Observations','Seasonal amplitude ratio', 'ENSO variance ratio')
axis(win), 
% make legend
h = [hp(Vanu) hp(Thierry) hp(Kim) hp(Woodroffe) hp(Helen) hp(Sandy) hp(Matthieu) hp(Robin)];
[legend_h,object_h,plot_h,text_strings] = columnlegend_h(2, h, lbl,'location','northeast','box','on');
set(plot_h(:),'LineStyle','none');
uistack(legend_h,'top'), set(legend_h,'position', [0.45   0.73    0.3800    0.1334])
 hold off
% PLOT MODELS
subplot(2,1,2), hold on
for m = 1:2*nm
    if m <= nm, mrkr{m} = 'v'; str = 'PI';
    else mrkr{m} = 's';    str = 'MH';
    end     
    for j = 1:3
        hm(m) = line(Xm(m,j),Ym(m,j),'MarkerFaceColor',colm(m,:),'marker',mrkr{m},'Markersize',8,'MarkerEdgeColor',bck);
        rectangle('Curvature', [1 1],'Position',[Xm(m,j)-Xm_iqr(m,j)/2 Ym(m,j)-Ym_iqr(m,j)/2 Xm_iqr(m,j) Ym_iqr(m,j)],'Edgecolor',colm(m,:),'Facecolor','none','LineWidth',1);
    end
end
axis(win)
%legends and axes
fancyplot_deco('General Circulation Models','Seasonal amplitude ratio', 'ENSO variance ratio')
[legend_h,object_h,plot_h,~] = columnlegend_h(3, hm(1:nm), models,'location','south','box','on');
set(plot_h(:),'LineStyle','none'); uistack(legend_h,'top'), set(legend_h,'position', [0.45   0.2    0.3800    0.1334])

% export
hepta_figprint('../figs/enso_var_vs_seas_ratios_BB_JEG_RAW',200);



% rescale variables by their expected errors
Xms = Xm./Xm_iqr; %Xms = Xms(:);
Yms = Ym./Ym_iqr; %Yms = Yms(:);
% fit TLS model
[Pm,Em,Xm_hats,Ym_hats] = fitTLS(Xms(:), Yms(:));
%Ym_hat = Ym_iqr(:).*Ym_hats; Xm_hat = Xm_iqr(:).*Xm_hats;

% boundaries
xm_min = 0; xm_max = 50;
ym_min = 1; ym_max = 3;

% boostrap estimate of regression coefficients
Pm_boot = bootstrp(Nb, @fitTLS,Xms(:),Yms(:));
Pm_se = std(Pm_boot(:,1));

xm = linspace(xm_min,xm_max,200); nxm = numel(xm);
Ym_hats = polyval(Pm,xm);
Lm = zeros(Nb,nxm); %matrix of regression lines
for b = 1:Nb
    Lm(b,:) = polyval(Pm_boot(b,:),xm);
end
ym = linspace(ym_min,ym_max,200); nym = numel(ym);
Lm_dens = zeros(nxm,nym);
for j = 1:nxm
    Lm_dens(j,:)= ksdensity(Lm(:,j),ym);
end
Lmq = quantile(Lm,[0.025 0.975],1);


% PLOT IT OUT
fig('ENSO_var vs seasonality'), clf
%set(gcf, 'Units','centimeters', 'Position',[0 0 9 7]) % 1 column width = 9cm
subplot(2,1,1)
% plot PDF for TLS fit
pcolor(xo,yo,Lo_dens'), shading interp
colormap(flipud(gray)); hold on, 
 % plot TLS fit
h2  = line(xo,Yo_hats,'color',bck,'linewidth',2);
hci  = plot(xo,Loq','k--'); w=1;
for s = 1:ns
    %hr(s)  = rectangle('Curvature', [1 1],'Position',[Xos(s)-w/2 Yos(s)-w/2 w w],'Edgecolor',col{s},'Facecolor','none','LineWidth',1);    
    hp(s)  = plot(Xos(s),Yos(s),'marker',mark(s),'MarkerFaceColor',col{s},'MarkerEdgeColor','k','MarkerSize',sz2(s));
end
fancyplot_deco('Proxy Observations','Scaled Seasonality', 'Scaled ENSO variance')
axis([0 20 -3 8]), 
labl{1} = ['TLS fit, \beta_o = ', sprintf('%+3.2g',Po(:,1)), '\pm',sprintf('%3.2g',2*Po_se)];
labl{2} = 'TLS fit, 95% CI';
hl = legend(gca,[h2 hci(1)],labl{:}), set(hl,style_l{:}), 
set(hl,'location','southeast','box','off'); hold off

% PLOT MODELS
subplot(2,1,2); hold on
% plot PDF for TLS fit
pcolor(xm,ym,Lm_dens'), shading interp
colormap(flipud(gray)); hold on,
 % plot TLS fit
h2  = line(xm,Ym_hats,'color',bck,'linewidth',2);
hci  = plot(xm,Lmq','k--');
for m = 1:2*nm
    if m <= nm, mrkr{m} = 'v'; str = 'PI';
    else mrkr{m} = 's';    str = 'MH';
    end     
    for j = 1:3
        hm(m) = line(Xms(m,j),Yms(m,j),'MarkerFaceColor',colm(m,:),'marker',mrkr{m},'Markersize',8,'MarkerEdgeColor',bck);
        %rectangle('Curvature', [1 1],'Position',[Xms(m,j)-w/2 Yms(m,j)-w/2 w w],'Edgecolor',colm(m,:),'Facecolor','none','LineWidth',1);
    end
end
axis([xm_min xm_max ym_min ym_max]); hold off
%legends and axes
for k = 1:nm
    lbl{k} = [models{k} ' PI'];
    lbl{nm+k} = [models{k} ' MH'];
end
fancyplot_deco('General Circulation Models','Scaled seasonality', 'Scaled ENSO variance')
% MAKE LEGENDS : requires duplicating axes
[legend_pi,object_pi,plot_pi,~] = columnlegend_h(3, hm(1:nm), lbl(1:nm),'location','northeast','box','off');
set(plot_pi(:),'LineStyle','none');
set(plot_h(:),'LineStyle','none'); uistack(legend_h,'top'), set(legend_h,'position', [0.45   0.2    0.3800    0.1334])

ah=axes('position',get(gca,'position'),'visible','off');
labl{1} = ['TLS fit, \beta_m = ', sprintf('%+3.2g',Pm(:,1)), '\pm',sprintf('%3.2g',2*Pm_se)];
labl{2} = 'TLS fit, 95% CI';
hl = legend(ah,[h2 hci(1)],labl{:}), set(hl,style_l{:}), 
set(hl,'location','southwest','box','off');
% export
hepta_figprint('../../../figs/Fig03',200);








% WAVELET ANALYSIS as in An & Choi 2013
% ==================
% wavelet parameters
Pad = 1; % pad the time series with zeroes (recommended)
dj  = 1/12; % this will do 12 sub-octaves per octave
S0  = 0.5;  % this says start at a scale of 0.5 years
Mother = 'Morlet';          
Cdelta = 0.776;   % scale factor for the MORLET wavelet
variance = 1;
continuous = find(~strcmp(cite,'Carre_2013'));

for r = continuous  % loop over records
    X = Ss(r).data; t = Ss(r).chron;  dt = Ss(r).res/12; n = length(t);
   %variance = var(X);
    
   [wave,period,scale,coi] = wavelet(X,dt,Pad,dj,S0,-1,Mother,-1);
   power = (abs(wave)).^2 ;      % compute wavelet power spectrum
   scaled_power = power ./ repmat(scale',[1 n]);   % [Eqn(24)]

   % Scale-average on relevant bands
   enso_band   = (scale >= 2 & scale <= 7);
   annual_band = (scale >= 0.8 & scale <=1.2);
   enso_power = variance*dj*dt/Cdelta*sum(scaled_power(enso_band,:));   % [Eqn(24)]
   annual_power = variance*dj*dt/Cdelta*sum(scaled_power(annual_band,:));   % [Eqn(24)]
   [rho_p(r),signif_p(r),pval_p(r)] = corr_sig(enso_power,annual_power);
   
   clab = ['\rho =  ', sprintf('%3.2g',rho_p(r)),', p = ',sprintf('%3.2g',pval_p(r))];
   
   ttl = ['Wavelet Power in ' Ss(r).meas ' at ' Ss(r).site,' (',Ss(r).reference,')'];
   fname = ['../figs/Rec_' sprintf('%02d',r) '_' Ss(r).sample_id '_enso_vs_AC.pdf'];
   % plot
%    fig('ENSO vs annual cycle'),clf
%    plot(t,enso_power,'b-',t,annual_power,'r-');
%    lbl{1}= 'ENSO band (2-7y)'; lbl{2}= 'annual band (0.8 - 1.2y)'; legend(lbl{:}), legend boxoff
%    fancyplot_deco(ttl,'Time(y)','Power')
%    text(0.1,0.9,clab,'Unit','Normalized')
%    export_fig(fname,'-r300','-cmyk')
end

   % pow = nextpow2(Sr(r).nyears);
   %  [wave,period,scale,coi] = wt_rect(TS, 'dj',dj,'MakeFigure','bw','AR1',g,'S0',2*dj,'MaxScale',2^pow);
   % [wave,period,scale,coi] = wt(TS, 'dj',dj,'MakeFigure','bw','AR1',g,'S0',2*dj,'MaxScale',2^pow);
   %    hcb = get(gcf, 'children');
   %    delete(hcb(1)), colormap(cejulien2),
   %    ttl = [Sr(r).site,' (',Sr(r).reference,') Morlet Wavelet Coefficients'];
   %    title(ttl,style_t{:},'FontName','Palatino');
   %    ylabel('Period (y)',style_l{:});
   %    set(gca,'XColor' , [.3 .3 .3], 'YColor', [.3 .3 .3], 'LineWidth', 1, 'TickDir','out');
   %    set(gca,'box', 'off'); xlabel('Time (y)')
   %    export_fig('../figs/Kiritimati_wavelet.pdf','-r300','-cmyk')


% ====================================
%  wavelet analysis of PMIP3 data
% ===================================
PMIP3 = load('../../../data/pmip3/coral_regional_avg.mat');

%exp = {'historical','past1000','midHolocene','piControl'}; 
exp = {'midHolocene','piControl'}; 

nexp = length(exp);
dt = 1/12.0; % all data are monthly

for ie = 1:nexp
    model_names{ie} = sort(fieldnames(eval(['PMIP3.' exp{ie}])));
    nm = length(model_names{ie});
    for im = 1:nm
        %extract timeseries
        struct = eval(['PMIP3.' exp{ie} '.' model_names{ie}{im}]);
        X = struct.center; n = length(X);
        t = [0:n-1]/12;
        % compute wavelet coefficients
        [wave,period,scale,coi] = wavelet(X,dt,Pad,dj,S0,-1,Mother,-1);
        power = (abs(wave)).^2 ;      % compute wavelet power spectrum
        scaled_power = power ./ repmat(scale',[1 n]);   % [Eqn(24)]
        variance = var(X);
        % Scale-average on relevant bands
        enso_band   = (scale >= 2 & scale <= 7);
        annual_band = (scale >= 0.8 & scale <=1.2);
        enso_power = variance*dj*dt/Cdelta*sum(scaled_power(enso_band,:));   % [Eqn(24)]
        annual_power = variance*dj*dt/Cdelta*sum(scaled_power(annual_band,:));   % [Eqn(24)]
        % correlation between ENSO and AC 
        [rho_m{ie}(im),signif_m{ie}(im),pval_m{ie}(im)] = corr_sig(enso_power,annual_power);       
        clab = ['\rho =  ', sprintf('%3.2g',rho_m{ie}(im)),', p = ',sprintf('%3.2g',pval_m{ie}(im))];
        ttl = [model_names{ie}{im} ' ' exp{ie}];
        fname = ['../../figs/wavelet/wavelet_' exp{ie} '_' model_names{ie}{im} '_center.pdf'];
        % plot
        fig('ENSO vs annual cycle'),clf
        subplot(2,1,1)
        plot(t,X);
        fancyplot_deco(['Pseudocoral expression of ' ttl],'Time(y)','Power')
        
        subplot(2,1,2)
        plot(t,enso_power,'b-',t,annual_power,'r-');
        lbl{1}= 'ENSO band (2-7y)'; lbl{2}= 'annual band (0.8 - 1.2y)'; legend(lbl{:}), legend boxoff
        fancyplot_deco([ 'Wavelet Power in ' ttl],'Time(y)','Power')
        text(0.1,0.9,clab,'Unit','Normalized')
        export_fig(fname,'-r200','-cmyk')
    end
end



clear X
% Plot correlations
fig('Correlation histograms'), clf
subplot(1,2,1)
sig_p = find(signif_p);
X.data = rho_p;
X.units = '\rho';
X.name = 'Proxy Obs';
[xh,H,kde]= KDE_hepta(X,10);
yl = get(gca,'YLim');
nsig = numel(sig_p);
% for j = 1:nsig
%     rh = rho_p(sig_p(j));
%     hl= plot([rh rh],yl,'k-.');
% end
% text(1.1*rh,1.4,' \leftarrow Palmyra modern','FontSize',12)
% the only record with a significant correlation is Palmyra modern, and
% this correlation is positive.

nm(1) = 0;
subplot(1,2,2)
col = hsv(4);
for ie = 1:nexp
    nm(ie+1) = length(rho_m{ie});
    rho_md(sum(nm(1:ie))+1:sum(nm(1:ie))+nm(ie+1)) = rho_m{ie};
    sig_md(sum(nm(1:ie))+1:sum(nm(1:ie))+nm(ie+1)) = signif_m{ie};
end
Y.data = rho_md;
Y.units = '\rho';
Y.name = 'PMIP3';
[xh,H,kde]= KDE_hepta(Y,10);
h = findobj(gca,'Type','patch');
set(h,'FaceColor',Lawn_green)
export_fig('../figs/AC_vs_ENSO_wavelet_rho_hist.pdf','-r200','-cmyk')
% only 13 out of 88 are significant. out of those , 6 are >0 and 7 are <0

  








