function scarpw=f_scarpwidth(xp,zp,npts,doplot,station)

%SCARPWIDTH Fault Scarp Width
%
% Description
%
%     Scarpwidth calculates and smooths the topographic slope along a profile 
%     using a moving mean window of specific size, and then fits a Gaussian
%     pulse function to the slope gradients. The fault scarp width is
%     estimated from the horizontal distance between the two peaks of the
%     double Gaussian function. This script is based in the method of:
%     Lu, et al., 2022, Frontiers in Earth Science, v. 10.
%
% Input arguments
%
%     xp,zp     topographic profile
%     npts      number of points of the moving mean window
%     doplot    plot results
%     station   site name
%
% Output arguments
%
%     scarpw    structure with results
%
% Author: Daniel Melnick (danmelski@gmail.com)
% Date: 2. April, 2025

%warning('off')

slp=diffxy(xp,zp); %slope
if npts>0
    sz=moving(slp,npts,'mean'); %moving mean
elseif npts==0
    sz=slp;
end
dslp=diffxy(xp,sz); %dz2/dx2 slope gradient

xxi=xp; xxi(isnan(dslp))=[];
dslp(isnan(dslp))=[];

sx=size(xxi); 
if sx(1)<sx(2)
    xxi=xxi';
end    

sy=size(dslp); 
if sy(1)<sy(2)
    dslp=dslp';
end    

tbl = table(xxi, dslp);
% Y = a + b*x + c*exp(-(x-d)^2/e) + d*exp(-(x-f)^2/g)

modelfun = @(b,x) b(1) + b(2) * x(:,1) + b(3) * exp(-(x(:,1) - b(4)).^2/b(5)) + b(6) * exp(-(x(:,1) - b(7)).^2/b(8));  
%initial guess
beta0 = [0, 0, 0, 10, 10, 0, 50, 20];   

mdl = fitnlm(tbl, modelfun, beta0);
coefficients = mdl.Coefficients{:, 'Estimate'};

X = linspace(min(xxi), max(xxi), 1000); 
Z = interp1(xp,zp,X,"linear");

yFitted = coefficients(1) + coefficients(2) * X + coefficients(3) * exp(-(X - coefficients(4)).^2 / coefficients(5)) + ...
	coefficients(6) * exp(-(X - coefficients(7)).^2 / coefficients(8));

miy=find(yFitted==min(yFitted));
may=find(yFitted==max(yFitted));

% Build output
scarpw.station=station;
scarpw.x=xp;
scarpw.z=zp;
scarpw.xxi=xxi;
scarpw.slope=slp;
scarpw.npts=npts;
scarpw.sslope=sz;
scarpw.dslope=dslp;
scarpw.X=X;
scarpw.Z=yFitted;
scarpw.x1=X(miy);
scarpw.x2=X(may);
scarpw.x1e=mdl.Coefficients.tStat(5);
scarpw.x2e=mdl.Coefficients.tStat(8);
scarpw.z1=Z(miy);
scarpw.z2=Z(may);
scarpw.ds1=yFitted(miy);
scarpw.ds2=yFitted(may);
scarpw.coef=mdl;
scarpw.rms=mdl.RMSE;
scarpw.r=mdl.Rsquared.Adjusted;
scarpw.w=X(miy)-X(may);
scarpw.we=sqrt(scarpw.x1e.^2+scarpw.x2e.^2);

% Plot
if doplot==1
mz=7;
figure(1), clf,
ax1=subplot(3,1,1);  hold on      
h1=plot(scarpw.x,scarpw.z,'.b');
errorbar(scarpw.x1,scarpw.z1,scarpw.x1e,'horizontal')
errorbar(scarpw.x2,scarpw.z2,scarpw.x2e,'horizontal')
h3=plot(scarpw.x1,scarpw.z1,'sr','markersize',mz,'MarkerFaceColor','r');
h2=plot(scarpw.x2,scarpw.z2,'dr','markersize',mz,'MarkerFaceColor','r');
legend([h1,h2,h3],'Puntos digitalizados','Base escarpe','Techo escarpe','Location','northwest')
ylabel('Altura (m)'); xlabel('Distancia horizontal (m)'); 
title(sprintf('%s Ancho del escarpe = %3.1f ± %2.1f m',scarpw.station,scarpw.w,scarpw.we)), box on
%
ax2=subplot(3,1,2); hold on
h1=plot(scarpw.x,scarpw.slope,'-k');
h2=plot(scarpw.x,scarpw.sslope,'-b');
legend([h1,h2],'Pendiente medida',sprintf('Media movil n=%u',npts),'Location','northwest')
ylabel('Pendiente (m/m)'); xlabel('Distancia horizontal (m)'); 
box on
%
ax3=subplot(3,1,3); hold on
h4=plot(scarpw.xxi,scarpw.dslope,'-r');
h1=plot(scarpw.X,scarpw.Z,'-b');
errorbar(scarpw.x1,scarpw.ds1,scarpw.x1e,'horizontal')
errorbar(scarpw.x2,scarpw.ds2,scarpw.x2e,'horizontal')
h2=plot(scarpw.x1,scarpw.ds1,'sr','markersize',mz,'MarkerFaceColor','r');
h3=plot(scarpw.x2,scarpw.ds2,'dr','markersize',mz,'MarkerFaceColor','r');
xlabel('Distancia horizontal (m)'); 
ylabel('Gradiente de pendiente')
legend([h4,h1,h2,h3],'d^2h/dx^2','Pulso Gaussiano','Minimo','Maximo','Location','southwest')

linkaxes([ax1 ax2 ax3], 'x'); box on
end

end
