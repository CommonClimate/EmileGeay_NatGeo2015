clear all
load ../data/proxy_db_holo_synthesis.mat

% FILTER based on errors
err = [Sr.enso_var_ratio_err]; 
% bad pupils err > 1 : Surprise Atoll, Tasmalum Sr/Ca and d18O
shortlist = find(err <= 1);  
%exclude = [52:55];  nr = length(Sr);
%shortlist = setdiff([1:nr],exclude);  

% PREPARE DATABASE TABLE
% ==================
Sr = Sr(shortlist);  nr = length(Sr);
age_mid = 1950-[Sr.age_mid]';
nyears = [Sr.nyears]';
res     = [Sr.res]';

lonp = round([Sr.lon]); lonp(lonp <0) = 360 + lonp(lonp <0);
latp = round([Sr.lat]);

archive ={Sr.archive}';
meas ={Sr.meas}';
site ={Sr.site}';
reso = [Sr.res]'; clear res
sample_id = {Sr.sample_id}';
modern_id = {Sr.modern_id}';
citekey  = {Sr.citekey};

for k=1:nr
   if latp(k) >=0
      slat{k} = ['$',num2str(latp(k)) '^{\circ}$N'];
   else
      slat{k} = ['$',num2str(-latp(k)) '^{\circ}$S'];
   end
   if (lonp(k) >=0 & lonp(k) <= 180)
      slon{k} = ['$',num2str(lonp(k)) '^{\circ}$E'];
   else
      slon{k} = ['$',num2str(round(360-lonp(k))) '^{\circ}$W'];
   end
   res{k} = int2str(reso(k));
   %
   range{k}= [int2str(Sr(k).year_i),'--',int2str(Sr(k).year_f)];
   if strcmp(citekey{k},'?')
      cite{k} = Sr(k).reference;
   else
      cite{k} = ['\citet{',citekey{k},'}'];
   end
   age{k} = int2str(age_mid(k));
   leng{k} = int2str(nyears(k));
   % TREAT ALL EXCEPTIONS INVOLVING _ or &
   sample_id{k} = strrep(sample_id{k},'_','\_');
   cite{k}      = strrep(cite{k},'&','\&');
   if strncmp('\delta^{18}O',meas{k},12)
      meas{k} = ['$' meas{k} '$'];
   end
end

rowLabels = cellstr(num2str([1:nr]'));
M=[archive site sample_id slat' slon' age' leng' res' meas cite'];
 

%  still a problem with Row Labels
cLabels = {'Record', 'Archive',' Site','Sample ID', 'Lat', 'Lon', 'Age', 'Length', 'Res', 'Meas','Reference'};
latextable(M,'Horiz',cLabels,'Vert',rowLabels,'name','../ancillary/proxy_db_table_input_raw.tex','Hline',[1]);
clear M  archive site sample_id slat slon age leng res meas cite range

% plot db stats
age_mid = [Sr.age_mid]';
nyears = [Sr.nyears]';
res     = [Sr.res]';
ligr = rgb('Silver');
fig('Database Stats'), clf
subplot(2,2,1)
x = [1:12]; hist(res,x);
h = findobj(gca,'Type','patch');
set(h,'FaceColor',ligr)
fancyplot_deco('a) Resolution','\Delta t (months)','Number of records',14);
subplot(2,2,2)
x = [10 25 35 50 75 100 125 150]'; hist(nyears,x)
hp = findobj(gca,'Type','patch'); set(hp,'FaceColor',ligr)
h = get(gca,'children');
set(h(1),'visible','off')
fancyplot_deco('b) Record length ','# of years','Number of records',14);
subplot(2,2,3:4)
x = [-8:1:2]*1000; xl = {'10','9' '8' '7' '6' '5' '4' '3' '2' '1' '0'}; %BP labels
hist(age_mid,x);  set(gca,'XTickLabel',xl(:))
h = findobj(gca,'Type','patch');
set(h,'FaceColor',ligr)
fancyplot_deco('c) Age coverage ','Median Age (ky BP)','Number of records',14);
export_fig('../figs/proxy_db_time_stats.pdf','-r200');

% Make the color figures for coral data anomaly plots, 6 to a page
tsPlot = 1;  % change this to 1 if you want to generate figures of all individual records
if tsPlot
   sites=char({Sr.site});
   %C
   for n=1:6:nr
      fig('Proxy Records'), clf
      %orient landscape
      %k=p;
      for k=n:n+5
         if k<=nr
            Sr(k).site
            subplot(3,2,k-n+1)
            tr=Sr(k).chron;  
            Xa=Sr(k).anom;
            dt = res(k);
            plot(tr-min(tr)+Sr(k).year_i,Xa,'r-'), grid; axis tight;
            ylabel(char(Sr(k).meas),'FontName','Palatino','FontSize',[10])
            xlabel('Time (years CE)','FontName','Palatino','FontSize',[10])
            % title
            rec=[num2str(k) ')'];
            obj=char(Sr(k).archive)
            meas=char(Sr(k).meas);
            lon_str=[num2str(round(Sr(k).lon)),'{}^{\circ}E'];
            lat3=Sr(k).lat;
            if (lat3>=0)
               lat_str=[num2str(round(lat3)),'{}^{\circ}N']
            else
               lat_str=[num2str(-round(lat3)),'{}^{\circ}S']
            end
            str=[rec,': ', obj,' ',meas,' at ', deblank(sites(k,:)),' (',lon_str,',',lat_str,'), \Delta t = ',int2str(dt)];
            title(str,'FontName','Palatino','FontSize',[10]); set(gca,'Box','Off','FontName','Palatino')
         end
      end
      filen=['../figs/proxy_database_', num2str(n),'_to_',num2str(n+5),'.pdf'];
      %pause
      export_fig(filen,'-cmyk','-r300');
   end
end

