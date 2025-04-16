clear 
load digitized_points.mat
station='PEN1';

% Loop window width
npts=2:2:14; numel(npts)
doplot=0; s=[];

for i=1:numel(npts)    
    scarpw=f_scarpwidth(geo_points(:,1),geo_points(:,2),npts(i),doplot,station);
    rm(i)=scarpw.rms;
    w(i)=scarpw.w;
    we(i)=scarpw.we;
end

% Plot relation
figure(1), clf, hold on
yyaxis left
plot(npts,rm,'o-b','MarkerFaceColor','auto')
title(sprintf('Sitio %s',station))
xlabel('Ancho ventana (n puntos)'), ylabel('RMSE')
yyaxis right
errorbar(npts,w,we)
plot(npts,w,'o-r','MarkerFaceColor','auto')
ylabel('Ancho del escarpe (m)'), box on
ylim([5 120])
xlim([1 15])
xw=15 ; yw=xw/1.7780; rect = [0,0,xw,yw]; set(gcf,'PaperUnits','centimeters','PaperType','A4','paperposition',rect);
mkdir('figs')
fout=sprintf('figs/%s_np_rmse.png',station); saveas(gcf,fout,'png');

% Plot 10 point analysis
npts=10;
scarpw=f_scarpwidth(geo_points(:,1),geo_points(:,2),npts,1,station);
xw=20 ; yw=25;
rect = [0,0,xw,yw]; set(gcf,'PaperUnits','centimeters','PaperType','A4','paperposition',rect);
fout=sprintf('figs/%s_prof_np%u.png',station,npts); saveas(gcf,fout,'png');
