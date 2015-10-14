clear all;
load ../../data/obs/proxy_db_holo_synthesis.mat
err = [So(:).enso_var_ratio_err];
export = 0
% define key arrays
ref = {So.reference};
site = {So.site};
archive = {So(:).archive}';
no = length(So);

c_lon = [So(:).lon]'; % fix longitudes:
c_lon(c_lon <0) = 360 + c_lon(c_lon <0);
c_lat = [So(:).lat]';
id_mod = {Sm.modern_id}; % find modern samples

% assign plotting parameters
for r = 1:no
    if strfind(site{r},'Vanuatu')
        So(r).col = rgb('SaddleBrown'); Vanuatu = r;
    elseif sum(strncmp(site{r},{'Baye','Surp'},4))>0
        So(r).col = rgb('Crimson');     Thierry = r;
    elseif strncmp(ref{r},{'Cobb'},4)
        So(r).col = rgb('LimeGreen');   Kim = r;
        %    elseif strfind(site{r},'Kiritimati')
        %       So(r).col = rgb('LimeGreen');   Kiri = r;
    elseif strncmp(ref{r},{'Wood'},4)
        So(r).col = rgb('LawnGreen');   Woodroffe = r;
    elseif strfind(ref{r},'Gregor')
        So(r).col  = rgb('DarkViolet');     Helen = r;
    elseif strfind(ref{r},'Tudhope')
        So(r).col  = rgb('DarkOrange'); Sandy = r;
    elseif strfind(ref{r},'Carr')
        So(r).col  = rgb('SteelBlue');   Matthieu = r;
    elseif strncmp(ref{r},'Welsh',5) | strncmp(ref{r},'Driscoll',8)
        So(r).col  = rgb('DeepSkyBlue'); Robin = r;
    end
    % assign plotting parameters
    if strcmp(archive{r},'coral')
        mark(r)='o';   sz(r) = 8;
    elseif strcmp(archive{r},'mollusk');
        mark(r)='p';   sz(r) = 12;
    end
end
col = {So.col};  sz2 = round(0.8*sz);
lbl{1} = 'Duprey et al. [2012]';
lbl{2} = 'this study';
lbl{3} = 'Cobb et al. [2013]';
lbl{4} = 'Woodroffe et al. [2003]';
lbl{5} = 'McGregor et al. [2004,2013]';
lbl{6} = 'Tudhope et al. [2001]';
lbl{7} = 'Carre et al. [2014]';
lbl{8} = 'Driscoll et al. [2014]';


% =================================================
%  SYNTHESIS FIGURES
% =================================================
load ../../data/obs/NCEP_OIv2_SODA_SSS_1981_2010_djf.mat
ilon = (lon>=100& lon <=300); ilat = (abs(lat)<=40);
lon_tp = lon(ilon); lat_tp = lat(ilat);
bck = rgb('Black');
style_l = {'FontName','Helvetica','Fontsize',8,'FontWeight','bold'};
% proxy locations + ENSO composites
% =================================================
fig('Proxy location'),clf
% set dimensions
set(gcf, 'Units','centimeters', 'Position',[0 0 9 7])
set(gcf, 'PaperPositionMode','auto')

m_proj('Robinson','clong',170,'lat',[-22 22],'lon',[130 300]);
hT = m_pcolor(lon_tp,lat_tp,pcoral_enso); caxis([-.4,.4]);
set(hT,'edgecolor','none'); hold on
colormap(cejulien2(21)); %colormap(pmkmp(17)); %
ci = linspace(-0.4,0.4,6);
m_coast('patch',rgb('Silver'));
m_grid('box','off','xtick',6,'ytick',7,'xlabeldir','middle', 'fontsize',5,'fontname','Helvetica');
hc = colorbar2('horiz','\delta^{18}O (permil)');
set(hc,'Position',[0.1300    0.250    0.7750    0.0407])
% add boxes.
h1 = m_line([120 120 180 180 120], [-20 0 0 -20 -20],'color',bck,'linewidth',1,'linestyle','-.');
h3 = m_line([270 270 280 280 270], [-10 0 0 -10 -10],'color',bck,'linewidth',1,'linestyle','-.');
h2 = m_line([190 190 240 240 190], [-5  5 5  -5  -5],'color',bck,'linewidth',1,'linestyle','-.');
% add proxies
for r=1:no
    hla(r)= m_line(c_lon(r),c_lat(r),'marker',mark(r),'MarkerEdgeColor',bck,'MarkerFaceColor',So(r).col,'Color',So(r).col,'linewidth',[1],'MarkerSize',sz(r),'linestyle','none');
    if strcmp(archive{r},'coral')
        h1 = hla(r); kc = r;
    elseif strcmp(archive{r},'mollusk');
        h2 = hla(r); km = r;
    end
end
% legend
h = [h1 h2];  lab = {'Coral', 'Mollusk'};
[legend_h,object_h,plot_h,text_strings] = columnlegend_h(2, h(:), lab,'location','SouthOutside','boxoff');
set(object_h(6),'MarkerFaceColor','w','linewidth',1),
set(object_h(4),'MarkerFaceColor','w','linewidth',1);
set(legend_h,style_l{:},'Position',[0.1300    0.25    0.5725    0.1237]),
ht = title(['Proxy locations & El Niño composite \delta^{18}O']);
set(ht,'FontName','Helvetica','Fontsize',10,'FontWeight','bold');

if export
    hepta_figprint('../../figs/Fig01',800);
    export_fig('../../figs/Fig01.eps','-cmyk','-r600','painters')
end

% ===================================
% INTERANNUAL AND SEASONAL VARIATIONS
% ===================================
seas = find(~isnan([So.seasonal_amp])); ns = numel(seas);
site_s = site(seas); ref_s  = ref(seas);

% total number of years covered
nyears_total = sum([So.nyears]);
age_mid = [So(:).age_mid];
year_i =  [So(:).year_i];
year_f =  [So(:).year_f];

% ========================================
% Stratification by longitude
% ========================================
reg_ind{1} = find(c_lon <=180);
reg_ind{2} = find(c_lon > 180 & c_lon < 240);
reg_ind{3} = find(c_lon > 240);
col = {So(:).col};
% load PMIP3 data
load ../../data/pmip3/PMIP3_ratios_distrib3.mat
% assign plotting parameters
nm = length(models);
%cmap = copper(nm);
%cmap = flipud(hsv(nm));
cmap = brewermap(nm,'Spectral');
%cmap = pmkmp(nm);
mdl = {'+','*','x','s','d','^','v','p','h'};
xl = {'9' '8' '7' '6' '5' '4' '3' '2' '1' '0'}; %BP labels
xt  = [-7:1:2]
rind = 1:20:numel(ratio_vec);
ylims = [0 3]; ymh = [0 3]; yls = [0 4]; ypi = [0 2.5];  yls_mh = [0 4];
yseas  =[0 3];
% plot the darn thing
fig('Regional breakdown'), clf
% set dimensions
set(gcf, 'Units','centimeters', 'Position',[0 0 17 12])
set(gcf, 'PaperPositionMode','auto')
p = panel();
p.pack([5 90 5]);
p(2).pack([1/3 1/3 1/3],[1/2 1/2]);
p.fontsize = 10; p.fontname = 'Helvetica';

for jr = 1:3
    % ENSO VARIANCE
    p(2,jr,1).pack('h',[1/6 1/6 2/3])
    p(2,jr,1,3).select()
    % proxies
    line([-9 2],[1 1],'linestyle','--','linewidth',1,'color',bck), hold on
    for r = reg_ind{jr}'
        loc = So(r).age_mid/1000;
        w   = So(r).nyears/1000;
        sig = So(r).age_std/1000;
        q   = So(r).enso_var_ratio_q; q(1) = max(q(1),0);
        line([So(r).year_i/1000 So(r).year_f/1000],[q(3) q(3)],'color',col{r},'linewidth',1)
        hr(r) = rectangle('Curvature', [1 1],'Position',[loc-2*sig q(1) 4*sig q(5)-q(1)],'Edgecolor',col{r},'Facecolor','none','LineWidth',1);
        hp(r) = plot(loc,q(3),'marker',mark(r),'MarkerFaceColor',col{r},'MarkerEdgeColor','k','MarkerSize',sz2(r));
    end
    hold off, axis([-8 2 ylims]), set(gca,'XTickLabel',xl(:),'XTick',xt);
    fancyplot_deco('','','',14)
    %% piControl GCMs
    p(2,jr,1,1).select(), hold on
    line([1 1],ypi,'linestyle','--','linewidth',1,'color',bck);
    for m = 1:nm
        f = sq(enso_var_ratio.kde.PIHT(m,jr,:));
        hpi_m(m) = plot(ratio_vec(rind),f(rind));
        set(hpi_m(m),'color',cmap(m,:),'linestyle','-','marker','none','linewidth',1);%mdl{m});
        [~,I] = min(abs(ratio_vec-enso_var_ratio.quant.PIHT(m,jr,3)));
        rmed = ratio_vec(I);
        line([rmed rmed],[0 f(I)],'color',cmap(m,:),'linestyle','-.','linewidth',1);
    end
    xlim(ylims); xlabel(region{jr},'fontsize',14,'fontweight','bold','fontname','Helvetica'),
    set(gca,'box','off'), ylim(ypi);
    set(gca,'Box', 'off', 'TickDir','out','TickLength',[.02 .02],'XMinorTick','on','YMinorTick','on', 'XGrid','on')
    set(gca,'XColor' , [.3 .3 .3], 'YColor', [.3 .3 .3], 'LineWidth', 1);
    view(-90,90), hold off
    
    
    %% MidHolocene GCMs
    p(2,jr,1,2).select(), hold on
    line([1 1],ymh,'linestyle','--','linewidth',1,'color',bck);
    for m = 1:nm
        f = sq(enso_var_ratio.kde.MHHT(m,jr,:));
        hmh_m(m) = plot(ratio_vec(rind),f(rind));
        set(hmh_m(m),'color',cmap(m,:),'linestyle','-','marker','none','linewidth',1);
        [~,I] = min(abs(ratio_vec-enso_var_ratio.quant.MHHT(m,jr,3)));
        rmed = ratio_vec(I);
        line([rmed rmed],[0 f(I)],'color',cmap(m,:),'linestyle','-.','linewidth',1);
    end
    view(90,-90), xlim(ylims), ylim(ymh);  hold off,
    set(gca,'XTickLabel',[]) % these are actually rotated Y ticks..
    set(gca,'Box', 'off', 'TickDir','out','TickLength',[.02 .02],'XMinorTick','on','YMinorTick','on', 'XGrid','on')
    set(gca,'XColor' , [.3 .3 .3], 'YColor', [.3 .3 .3], 'LineWidth', 1);
    
    % SEASONALITY
    p(2,jr,2).pack('h',[1/6 1/6 2/3])
    p(2,jr,2,3).select()
    %proxies
    hl = line([-8 2],[1 1],'linestyle','--','linewidth',1,'color',bck), hold on
    seas_reg = intersect(reg_ind{jr},seas)';
    for r = seas_reg
        loc = So(r).age_mid/1000;
        w   = So(r).nyears/1000;
        sig = So(r).age_std/1000;
        q   = So(r).seasonal_amp_ratio_q;
        line([So(r).year_i/1000 So(r).year_f/1000],[q(3) q(3)],'color',col{r},'linewidth',1)
        hr(r) = rectangle('Curvature', [1 1],'Position',[loc-2*sig q(1) 4*sig q(5)-q(1)],'Edgecolor',col{r},'Facecolor','none','LineWidth',1);
        hp(r) = plot(So(r).age_mid/1000,q(3),'marker',mark(r),'MarkerFaceColor',col{r},'MarkerEdgeColor','k','MarkerSize',sz2(r));
    end
    hold off, axis([-8 2 yseas(1) yseas(2)]); set(gca,'XTickLabel',xl(:),'XTick',xt)
    fancyplot_deco('','','',10)
    %% piControl GCMs
    p(2,jr,2,1).select(), hold on
    line([1 1],yls,'linestyle','--','linewidth',1,'color',bck);
    for m = 1:nm
        f = sq(seas_amp_ratio.kde.PIHT(m,jr,:));
        hpi_m(m) = plot(ratio_vec(rind),f(rind));
        set(hpi_m(m),'color',cmap(m,:),'linestyle','-','marker','none','linewidth',1);
        [~,I] = min(abs(ratio_vec-seas_amp_ratio.quant.PIHT(m,jr,3)));
        rmed = ratio_vec(I);
        line([rmed rmed],[0 f(I)],'color',cmap(m,:),'linestyle','-.','linewidth',1);
    end
    xlim(yseas); %xlabel(region{jr},'fontsize',14,'fontweight','bold'),
    set(gca,'box','off'),
    view(-90,90),ylim(yls);  hold off
    set(gca,'Box', 'off', 'TickDir','out','TickLength',[.02 .02],'XMinorTick','on','YMinorTick','on', 'XGrid','on')
    set(gca,'XColor' , [.3 .3 .3], 'YColor', [.3 .3 .3], 'LineWidth', 1);
    
    %% MidHolocene GCMs
    p(2,jr,2,2).select(), hold on
    line([1 1],yls_mh,'linestyle','--','linewidth',1,'color',bck),
    for m = 1:nm
        f = sq(seas_amp_ratio.kde.MHHT(m,jr,:));
        hmh_m(m) = plot(ratio_vec(rind),f(rind));
        set(hmh_m(m),'color',cmap(m,:),'linestyle','-','marker','none','linewidth',1);
        [~,I] = min(abs(ratio_vec-seas_amp_ratio.quant.MHHT(m,jr,3)));
        rmed = ratio_vec(I);
        line([rmed rmed],[0 f(I)],'color',cmap(m,:),'linestyle','-.','linewidth',1);
    end
    xlim(yseas); view(90,-90), ylim(yls_mh);   hold off, title('');
    set(gca,'XTickLabel',[]) % these are actually rotated Y ticks..
    set(gca,'Box', 'off', 'TickDir','out','TickLength',[.02 .02],'XMinorTick','on','YMinorTick','on', 'XGrid','on')
    set(gca,'XColor' , [.3 .3 .3], 'YColor', [.3 .3 .3], 'LineWidth', 1);
end
%  FIX LABELS
p(2,1,1,1).select()
title('PI','fontsize',12, 'fontweight','bold');
p(2,1,1,2).select()
title('MH','fontsize',12, 'fontweight','bold');
p(2,1,1,3).select()
title('Observations','fontsize',12, 'fontweight','bold');
% main titles
text(-10,4,'ENSO variance ratio','fontsize',14, 'fontweight','bold');
text(10,4,'AC amplitude ratio','fontsize',14, 'fontweight','bold');
%
p(2,1,2,1).select()
title('PI','fontsize',12, 'fontweight','bold');
p(2,1,2,2).select()
title('MH','fontsize',12, 'fontweight','bold');
p(2,1,2,3).select()
title('Observations','fontsize',12, 'fontweight','bold');
%
p(2,3,1,1).select()
ylabel('Density','fontsize',10);
p(2,3,1,3).select()
fancyplot_deco('','Time (ky BP)','',10);
p(2,3,1,2).select()
ylabel('Density');
p(2,3,2,1).select()
ylabel('Density');
p(2,3,2,3).select()
fancyplot_deco('','Time (ky BP)','',10);
p(2,3,2,2).select()
ylabel('Density');

% make legend
p(3).select()
set(gca,'visible','off')
[legend_h,object_h,plot_h,text_strings] = columnlegend_h(4, hpi_m, models,'location','Best','box','on');
uistack(legend_h,'top')
set(object_h(9:24),'LineWidth',2);
set(legend_h,'position', [0.15   -0.06    0.3800    0.1334])

% print out
if export
    hepta_figprint('../../figs/Fig02',800);
end
% cube_helix colormap
%cmap = cubehelix(10,.5,-1.5,3,1,[.2,1],[0,1])
% fig, hold on
% for m = 1:nm
%     f = sq(enso_var_ratio.kde.PIHT(m,jr,:));
%     hpi_m(m) = plot(ratio_vec(rind),f(rind));
%     set(hpi_m(m),'color',cmap(m,:),'linestyle','-','marker','none','linewidth',2);%mdl{m});
%     [~,I] = min(abs(ratio_vec-enso_var_ratio.quant.PIHT(m,jr,3)));
%     rmed = ratio_vec(I);
%     line([rmed rmed],[0 f(I)],'color',cmap(m,:),'linestyle','-.','linewidth',1);
% end
% xlim(ylims);
% set(gca,'box','off'), ylim(ylm);
% view(-90,90), hold off
% [legend_h,object_h,plot_h,text_strings] = columnlegend_h(4, hpi_m, models,'location','Best','box','on');
% uistack(legend_h,'top')
% set(object_h(9:24),'LineWidth',2);
% hepta_figprint('../../figs/colorscheme_prism',200);



% Probability exercise:
% ====================
% central Pacific reduction
no = length(So);
for jr = 1:3
    nreg = numel(reg_ind{jr});
    loc = zeros(nreg,1); var_ratio = zeros(nreg,1000);
    for r = 1:nreg
        loc(r) = So(reg_ind{jr}(r)).age_mid/1000;
        var_ratio(r,:) = So(reg_ind{jr}(r)).enso_var_ratio_b;
        %weight(r)    = 1/So(r).enso_var_ratio_err;
    end
    idx = (loc >= -3 & loc <= -1);
    v = var_ratio(idx,:);
    enso_var_ratio_3k5k(jr,:) =  quantile(v(:),[0.025 0.5 0.975]);
    enso_var_red_3k5k(jr,:) = 100*(quantile(1-v(:),[0.025 0.5 0.975]));
    %
    idx = (loc >= -5.5 & loc <= -3.5);  v = var_ratio(idx,:);
    enso_var_ratio_5k7k(jr,:) =  quantile(v(:),[0.025 0.5 0.975]);
    enso_var_red_5k7k(jr,:) = 100*(quantile(1-v(:),[0.025 0.5 0.975]));
    
    %enso_var_ratio_5k7k(jr) = sum(var_ratio(idx).*weight(idx))/sum(weight(idx));
end
% E pacific (Peruvian bivalves)
bivalves_age = [So(reg_ind{3}).age_mid]/1000;

enso_var_red_3k5k(3,:) = sort(100*(1-So(reg_ind{jr}(3)).enso_var_ratio_q([1,3,5])));
enso_var_red_5k7k(3,:) = sort(100*(1-So(reg_ind{jr}(4)).enso_var_ratio_q([1,3,5])));
%
for j = 1:3
    for k = 1:3
        red_table_5k7k{j,k} =  [num2str(enso_var_red_5k7k(j,k),'%3.0f'), '\%'];
        red_table_3k5k{j,k} =  [num2str(enso_var_red_3k5k(j,k),'%3.0f'), '\%'];
    end
end

latextable(red_table_5k7k,'name','red_table_5k7k.tex','Hline',1,'Vline',1);
latextable(red_table_3k5k,'name','red_table_3k5k.tex','Hline',1,'Vline',1);



for m = 1:nm
    % under a PI null
    f = sq(enso_var_ratio.kde.PIHT(m,2,:));
    idx = (ratio_vec <= enso_var_ratio_3k5k(2,2));
    table{1,m} = [num2str(100*trapz(ratio_vec(idx),f(idx)),'%3.2f'),'\%'];
    % under an MH null
    f = sq(enso_var_ratio.kde.MHHT(m,2,:));
    table{2,m} = [num2str(100*trapz(ratio_vec(idx),f(idx)),'%3.2f'),'\%'];
end

cLab = {'$H_0$', models{:}}
rLab = {'PI null','MH null'};
latextable(table,'Horiz',cLab','Vert',rLab','name','probability_observing_3k5k.tex','Hline',1,'Vline',1);


save ../data/proxy_db_holo_synthesis.mat cmap lbl -append


% plot the darn thing
% fig('Regional breakdown'), clf
% p = panel();
% p.pack([4/5 1/5]);
% p(1).pack([1/2 1/2],[1/3 1/3 1/3]);
% p.fontsize = 8; p.fontname = 'Helvetica';
% %p.select('all');
% %p.identify()
%
% for jr = 1:3
%     % ENSO VARIANCE
%     p(1,1,jr).pack('h',[1/6 2/3 1/6])
%     p(1,1,jr,2).select()
%     % proxies
%     line([-9 2],[1 1],'linestyle','--','linewidth',1,'color',bck), hold on
%     for r = reg_ind{jr}'
%         loc = So(r).age_mid/1000;
%         w   = max(So(r).age_std,So(r).nyears/2)/1000;
%         q   = So(r).enso_var_ratio_q;
%         hr(r) = rectangle('Curvature', [1 1],'Position',[loc-2*w q(1) 4*w q(5)-q(1)],'Edgecolor',col{r},'Facecolor','none','LineWidth',1);
%         hp(r) = plot(So(r).age_mid/1000,So(r).enso_var_ratio,'marker',mark(r),'MarkerFaceColor',col{r},'MarkerEdgeColor','k','MarkerSize',sz2(r));
%     end
%     hold off, axis([-9 2 ylims]), set(gca,'XTickLabel',xl(:))
%     %ttl = ['ENSO variance, ' region{jr}];
%     fancyplot_deco(region{jr},'Time (ky BP)','',10)
%     %% piControl GCMs
%     p(1,1,jr,1).select(), hold on
%     for m = 1:nm
%         f = sq(kde.PIHT(m,jr,:));
%         hpi_m(m) = plot(ratio_vec(rind),f(rind));
%         set(hpi_m(m),'color',cmap(m,:),'linestyle','-','marker',mdl{m});
%         [~,I] = min(abs(ratio_vec-var_ratio_q.PIHT(m,jr,3)));
%         rmed = ratio_vec(I);
%         line([rmed rmed],[0 f(I)],'color',cmap(m,:),'linestyle','-','linewidth',2);
%     end
%     xlim(ylims); xlabel('ENSO variance Ratio','fontsize',12), ylabel('Density');
%     title('piControl');
%     view(-90,90), hold off
%
%     %% MidHolocene GCMs
%     p(1,1,jr,3).select(), hold on
%     for m = 1:nm
%         f = sq(kde.MHHT(m,jr,:));
%         hmh_m(m) = plot(ratio_vec(rind),f(rind));
%         set(hmh_m(m),'color',cmap(m,:),'linestyle','-','marker',mdl{m});
%         [~,I] = min(abs(ratio_vec-var_ratio_q.MHHT(m,jr,3)));
%         rmed = ratio_vec(I);
%         line([rmed rmed],[0 f(I)],'color',cmap(m,:),'linestyle','-','linewidth',2);
%     end
%     view(90,-90), hold off, title('midHolocene'); ylabel('Density');
%
%     % SEASONALITY
%     p(1,2,jr).pack('h',[1/6 2/3 1/6])
%     p(1,2,jr,2).select()
%     %proxies
%     hl = line([-9 2],[1 1],'linestyle','--','linewidth',1,'color',bck), hold on
%     seas_reg = intersect(reg_ind{jr},seas);
%     for r = seas_reg
%         loc = So(r).age_mid/1000;
%         w   = max(So(r).age_std,So(r).nyears/2)/1000;
%         q   = So(r).seasonal_amp_ratio_q;
%         hr(r) = rectangle('Curvature', [1 1],'Position',[loc-2*w q(1) 4*w q(5)-q(1)],'Edgecolor',col{r},'Facecolor','none','LineWidth',1);
%         hp(r) = plot(So(r).age_mid/1000,So(r).seasonal_amp_ratio,'marker',mark(r),'MarkerFaceColor',col{r},'MarkerEdgeColor','k','MarkerSize',sz2(r));
%     end
%     hold off, axis([-9 2 0 4]); set(gca,'XTickLabel',xl(:))
%     fancyplot_deco(region{jr},'Time (ky BP)','',10)
%     %% piControl GCMs
%     p(1,2,jr,1).select(), hold on
%     for m = 1:nm
%         f = sq(kde.PIHT(m,jr,:));
%         hpi_m(m) = plot(ratio_vec(rind),f(rind));
%         set(hpi_m(m),'color',cmap(m,:),'linestyle','-','marker',mdl{m});
%         [~,I] = min(abs(ratio_vec-var_ratio_q.PIHT(m,jr,3)));
%         rmed = ratio_vec(I);
%         line([rmed rmed],[0 f(I)],'color',cmap(m,:),'linestyle','-','linewidth',2);
%     end
%     xlim(ylims); xlabel('Seasonal Amplitude','fontsize',12), ylabel('Density');
%     title('piControl');
%     view(-90,90), hold off
%
%     %% MidHolocene GCMs
%     p(1,2,jr,3).select(), hold on
%     for m = 1:nm
%         f = sq(kde.MHHT(m,jr,:));
%         hmh_m(m) = plot(ratio_vec(rind),f(rind));
%         set(hmh_m(m),'color',cmap(m,:),'linestyle','-','marker',mdl{m});
%         [~,I] = min(abs(ratio_vec-var_ratio_q.MHHT(m,jr,3)));
%         rmed = ratio_vec(I);
%         line([rmed rmed],[0 f(I)],'color',cmap(m,:),'linestyle','-','linewidth',2);
%     end
%     view(90,-90), hold off, title('midHolocene'); ylabel('Density');
% end
%
% hepta_figprint('../../figs/variance_ratios_regional_kde',400);







% OLD WAY
% fig('Regional breakdown'), clf
% xl = {'10' '8' '6' '4' '2' '0'}; %BP labels
% for jr = 1:3
%     % ENSO VARIANCE
%     subplot(2,3,jr)
%     line([-9 2],[1 1],'linestyle','--','linewidth',1,'color',bck), hold on
%     for m = 1:nm
%         % add PMIP3 block-boostrap quantiles levels as lines
%         % PI/HT levels
%         line([-9 2],[var_ratio_q.PIHT(m,jr,1) var_ratio_q.PIHT(m,jr,1)],'linestyle','-.','linewidth',1,'color',cmap(m,:));
%         line([-9 2],[var_ratio_q.PIHT(m,jr,3) var_ratio_q.PIHT(m,jr,3)],'Marker','s','linewidth',2,'color',cmap(m,:));
%         line([-9 2],[var_ratio_q.PIHT(m,jr,5) var_ratio_q.PIHT(m,jr,5)],'linestyle','-','linewidth',1,'color',cmap(m,:));
%         % MH/HT level
%         line([-7 -4],[var_ratio_q.MHHT(m,jr,1) var_ratio_q.MHHT(m,jr,1)],'linestyle','-.','linewidth',1,'color',cmap(m,:));
%         line([-7 -4],[var_ratio_q.MHHT(m,jr,3) var_ratio_q.MHHT(m,jr,3)],'Marker','s','linewidth',2,'color',cmap(m,:));
%         line([-7 -4],[var_ratio_q.MHHT(m,jr,5) var_ratio_q.MHHT(m,jr,5)],'linestyle','-','linewidth',1,'color',cmap(m,:));
%     end
%
%     for r = reg_ind{jr}'
%         loc = So(r).age_mid/1000;
%         w   = max(So(r).age_std,So(r).nyears/2)/1000;
%         q   = So(r).enso_var_ratio_q;
%         hr(r) = rectangle('Curvature', [1 1],'Position',[loc-2*w q(1) 4*w q(5)-q(1)],'Edgecolor',col{r},'Facecolor','none','LineWidth',1);
%         hp(r) = plot(So(r).age_mid/1000,So(r).enso_var_ratio,'marker',mark(r),'MarkerFaceColor',col{r},'MarkerEdgeColor','k','MarkerSize',sz2(r));
%     end
%     hold off, axis([-9 2 0 4]), set(gca,'XTickLabel',xl(:))
%     ttl = ['ENSO variance, ' region{jr}]; fancyplot_deco(ttl,'Time (ky BP)','Ratio',14)
%
%     % SEASONALITY
%     subplot(2,3,jr+3)
%     hl = line([-9 2],[1 1],'linestyle','--','linewidth',1,'color',bck), hold on
%     for m = 1:nm
%         % add PMIP3 block-boostrap quantiles levels as lines
%         % PI/HT levels
%         hpi_l(m) = line([-9 2],[var_ratio_q.PIHT(m,jr,1) var_ratio_q.PIHT(m,jr,1)],'linestyle','-.','linewidth',1,'color',cmap(m,:));
%         hpi_m(m) = line([-9 2],[var_ratio_q.PIHT(m,jr,3) var_ratio_q.PIHT(m,jr,3)],'Marker','s','linewidth',2,'color',cmap(m,:));
%         hpi_u(m) = line([-9 2],[var_ratio_q.PIHT(m,jr,5) var_ratio_q.PIHT(m,jr,5)],'linestyle','-','linewidth',1,'color',cmap(m,:));
%         lpi{m} = [models{m} ', piControl'];
%         % MH/HT level
%         hmh_l(m) = line([-7 -4],[var_ratio_q.MHHT(m,jr,1) var_ratio_q.MHHT(m,jr,1)],'linestyle','-.','linewidth',1,'color',cmap(m,:));
%         hmh_m(m) = line([-7 -4],[var_ratio_q.MHHT(m,jr,3) var_ratio_q.MHHT(m,jr,3)],'Marker','s','linewidth',2,'color',cmap(m,:));
%         hmh_u(m) = line([-7 -4],[var_ratio_q.MHHT(m,jr,5) var_ratio_q.MHHT(m,jr,5)],'linestyle','-','linewidth',1,'color',cmap(m,:));
%         lmh{m} = [models{m} ', Mid-Holocene'];
%     end
%     seas_reg = intersect(reg_ind{jr},seas);
%     for r = seas_reg
%         loc = So(r).age_mid/1000;
%         w   = max(So(r).age_std,So(r).nyears/2)/1000;
%         q   = So(r).seasonal_amp_ratio_q;
%         hr(r) = rectangle('Curvature', [1 1],'Position',[loc-2*w q(1) 4*w q(5)-q(1)],'Edgecolor',col{r},'Facecolor','none','LineWidth',1);
%         hp(r) = plot(So(r).age_mid/1000,So(r).seasonal_amp_ratio,'marker',mark(r),'MarkerFaceColor',col{r},'MarkerEdgeColor','k','MarkerSize',sz2(r));
%     end
%
%     hold off, axis([-9 2 0 4]); set(gca,'XTickLabel',xl(:))
%     ttl = ['Seasonality, ' region{jr}];
%     fancyplot_deco(ttl,'Time (ky BP)','Ratio',14);
% end
%
% subplot(2,3,6)
% h = [hpi_m hmh_m];  lab = horzcat(lpi, lmh);
% [legend_h,object_h,plot_h,text_strings] = columnlegend_h(2, h, lab,'location','Best','box','on');
% set([hpi_l' hpi_m' hpi_u' hmh_l' hmh_m' hmh_u'] ,'Visible','off'), set(hl,'Visible','off'); % erase lines
% set(gca,'visible','off')
%
% hepta_figprint('../../figs/variance_ratios_regional_quantiles');
%









% ====================================
% Can any millennium be singled out ?
% ====================================

% boostrap over samples to obtain the distribution


% ============================================================
%   Link to insolation
% ============================================================
% Generate insolation field using Huybers & Eisenman method
kyr = 0:0.2:21;
day = 1:5:365;
lat = 0;

[kyr3 day3 lat3] = ndgrid(kyr,day,lat);
[ny, nd ,nl]  = size(kyr3);
Fsw           = daily_insolation(kyr3,lat3,day3);
%Fsws = Fsw.*cosd(lat);   may be unnecessary since Earth's curvature is
%already taken into account.
%
% % reshape
% Fsw2d = reshape(Fsw,[ny nd*nl]);
% Z = center(Fsw2d);
% % conduct SVD
% [P,S,E]=svd(Z);
% q = min(size(S));
% % rewrite the large eigenvalue matrix to a vector and
% % apply the appropriate normalization
% eigv=diag(S(1:q,1:q)).^2/(q-1);    % has q elements
% spec=100*diag(eigv)/sum(diag(eigv)); % percentage of variance explained
%
% % Compute Sept-March
% Sep = squeeze(sum(Fsw(datenum(0,9,1)-1:datenum(0,9,30)-1,:,:)));
% March = squeeze(sum(Fsw(datenum(0,3,1)-1:datenum(0,3,31)-1,:,:)));
% % --->  some very weird mid-Holocene glitches, either in EOFs or direct
% % computations.










% % investigate outliers
% R = [So.enso_var_ratio];
% fig('pathology'), clf
% rp = 32; id = 10;
% to = So(rp).chron - min(So(rp).chron);
% tm = Sm(id).chron - min(Sm(id).chron);
% plot(to,So(rp).anom,'r-',tm,Sm(id).anom,'k-')

% 1) Woodruffe 2003 CW-2 sample : 23 y
% 2) Correge Surprise SU3 (17y) and SUmod (8) ---> too short
% 3) Correge Tasmalum  31y  vs KIL04 18O 23y


