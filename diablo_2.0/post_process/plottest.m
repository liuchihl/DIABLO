figure;
 set(gcf,'units','centimeters','paperunits','centimeters')
 set(gcf,'PaperType','A4');
 pp=[0.63 .9 26 24];
ps=[0 0 pp(3)/1.1 pp(4)/1.1];
 set(gcf,'paperposition',pp);
 set(gcf,'position',ps);

 i=1;
 ax1 = axes('position',[.08 .725 .26 .22]);
   [cc,hh] = contourf(x{i},y{i},squeeze(K3d{i})',100);
set(hh,'edgecolor','none');
title(' $K_{3d}$','fontname','times','fontsize',18,...
    'fontweight','normal','interpreter','latex');
%caxis([0 0.04]);
% caxis([-1 1]);
colormap(ax1,'default')
text(1.5,-2.5,'$d=10, t=188$','interpreter','latex',...
'fontsize',15,'color','w');
ylabel('Z','fontsize',12);
ylim([y{i}(1)+10-3.29 y{i}(1)+10+3.29]);
set(gca,'fontsize',12,'tickdir','out',...
    'ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;


ax2 = axes('position',[.35 .725 .26 .22]);
  [cc,hh] = contourf(x{i},y{i},squeeze(SP_background{i})',100);
set(hh,'edgecolor','none');
%caxis([-.002 .01]);
colormap(ax2,jet);
title('$SP_{background}$','interpreter','latex',...
'fontsize',18,'color','k');
set(gca,'fontsize',12,'yticklabel',[],'tickdir','out',...
    'ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;
ylim([y{i}(1)+10-3.29 y{i}(1)+10+3.29]);


ax3 = axes('position',[.62 .725 .26 .22]);
  [cc,hh] = contourf(x{i},y{i},squeeze(D3d{i})',100);set(hh,'edgecolor','none');
colormap(ax3,flipud((cbrewer('seq', 'YlGn', 100))));
caxis([-3.5,0]*1e-3);
title('$\epsilon_{3d}$','interpreter','latex',...
'fontsize',15,'color','k');
set(gca,'fontsize',12,'tickdir','out','ticklength',[0.02 0.02],...
    'xminortick','on','yminortick','on');box on;
ylim([y{i}(1)+10-3.29 y{i}(1)+10+3.29]);
set(gca,'fontsize',12,'yticklabel',[],'tickdir','out',...
    'ticklength',[0.02 0.02],'xminortick','on','yminortick','on');



i=2;
ax1 = axes('position',[.08 .45 .26 .22]);
[cc,hh] = contourf(x{i},y{i},squeeze(K3d{i})',100);
set(hh,'edgecolor','none');
%caxis([0 0.04]);
% caxis([-1 1]);
colormap(ax1,'default')
text(1.5,-2.5,'$d=10, t=188$','interpreter','latex',...
'fontsize',15,'color','w');
ylabel('Z','fontsize',12);
ylim([y{i}(1)+10-3.29 y{i}(1)+10+3.29]);
set(gca,'fontsize',12,'tickdir','out',...
    'ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;


ax2 = axes('position',[.35 .45 .26 .22]);
    [cc,hh] = contourf(x{i},y{i},squeeze(SP_background{i})',100);
set(hh,'edgecolor','none');
%caxis([-.002 .01]);
colormap(ax2,jet);
set(gca,'fontsize',12,'yticklabel',[],'tickdir','out',...
    'ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;
ylim([y{i}(1)+10-3.29 y{i}(1)+10+3.29]);


ax3 = axes('position',[.62 .45 .26 .22]);
  [cc,hh] = contourf(x{i},y{i},squeeze(D3d{i})',100);set(hh,'edgecolor','none');
colormap(ax3,flipud((cbrewer('seq', 'YlGn', 100))));
caxis([-3.5,0]*1e-3);
set(gca,'fontsize',12,'tickdir','out','ticklength',[0.02 0.02],...
    'xminortick','on','yminortick','on');box on;
ylim([y{i}(1)+10-3.29 y{i}(1)+10+3.29]);
set(gca,'fontsize',12,'yticklabel',[],'tickdir','out',...
    'ticklength',[0.02 0.02],'xminortick','on','yminortick','on');

i=3;
 ax1 = axes('position',[.08 .13 .26 .22]);
[cc,hh] = contourf(x{i},y{i},squeeze(K3d{i}'),100);
set(hh,'edgecolor','none');
%caxis([0 0.03]);
colormap(ax1,'default')
ylabel('Z','fontsize',12);xlabel('X','fontsize',12);
col = colorbar('location','southoutside','position',[.09,0.056 0.19 0.015]);
set(gca,'fontsize',12,...
    'tickdir','out','ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;
text(1.5,-5.8,'$d=3, t=266$','interpreter','latex',...
'fontsize',15,'color','w');
ylim([y{i}(1) 0])
 
ax2 = axes('position',[.35 .13 .26 .22]);
[cc,hh] = contourf(x{i},y{i},squeeze(SP_background{i})',100);
set(hh,'edgecolor','none');
%caxis([-4e-3 0.01]);
colormap(ax2,jet);
xlabel('X','fontsize',12);
set(gca,'fontsize',12,'yticklabel',[],...
    'tickdir','out','ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;
ylim([y{i}(1) 0]);


ax3 = axes('position',[.62 .13 .26 .22]);
[cc,hh] = contourf(x{i},y{i},squeeze(D3d{i})',100);set(hh,'edgecolor','none');
% hold on; [cc,hh]=contour(x{i},y{i},tke_3D_m',10,'color','k');%set(cc,'color','k');
colormap(ax3,flipud((cbrewer('seq', 'YlGn', 100))));
col = colorbar('location','southoutside','position',[.75,0.056 0.19 0.015]);
xlabel('X','fontsize',12);
set(gca,'fontsize',12,'yticklabel',[],...
    'tickdir','out','ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
ylim([y{i}(1) 0]);

%% plot SP_2d, SP_sheardeformation and BF3d
figure;
 set(gcf,'units','centimeters','paperunits','centimeters')
 set(gcf,'PaperType','A4');
 pp=[0.63 .9 26 24];
ps=[0 0 pp(3)/1.1 pp(4)/1.1];
 set(gcf,'paperposition',pp);
 set(gcf,'position',ps);

 i=1;
ax1 = axes('position',[.08 .725 .26 .22]);
[cc,hh] = contourf(x{i},y{i},squeeze(SP_2d{i})',100);
set(hh,'edgecolor','none');
title(' $SP_{KH}$','fontname','times','fontsize',18,...
    'fontweight','normal','interpreter','latex');
%caxis([0 0.04]);
% caxis([-1 1]);
colormap(ax1,'default')
text(1.5,-2.5,'$d=10, t=188$','interpreter','latex',...
'fontsize',15,'color','w');
ylabel('Z','fontsize',12);
ylim([y{i}(1)+10-3.29 y{i}(1)+10+3.29]);
set(gca,'fontsize',12,'tickdir','out',...
    'ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;


ax2 = axes('position',[.35 .725 .26 .22]);
  [cc,hh] = contourf(x{i},y{i},squeeze(SP_sheardeform{i})',100);
set(hh,'edgecolor','none');
%caxis([-.002 .01]);
colormap(ax2,jet);
title('$SP_{shear}$','interpreter','latex',...
'fontsize',18,'color','k');
set(gca,'fontsize',12,'yticklabel',[],'tickdir','out',...
    'ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;
ylim([y{i}(1)+10-3.29 y{i}(1)+10+3.29]);


ax3 = axes('position',[.62 .725 .26 .22]);
  [cc,hh] = contourf(x{i},y{i},squeeze(BF3d{i})',100);set(hh,'edgecolor','none');
colormap(ax3,flipud((cbrewer('seq', 'YlGn', 100))));
caxis([-3.5,0]*1e-3);
title('$BF_{3d}$','interpreter','latex',...
'fontsize',15,'color','k');
set(gca,'fontsize',12,'tickdir','out','ticklength',[0.02 0.02],...
    'xminortick','on','yminortick','on');box on;
ylim([y{i}(1)+10-3.29 y{i}(1)+10+3.29]);
set(gca,'fontsize',12,'yticklabel',[],'tickdir','out',...
    'ticklength',[0.02 0.02],'xminortick','on','yminortick','on');



i=2;
ax1 = axes('position',[.08 .45 .26 .22]);
[cc,hh] = contourf(x{i},y{i},squeeze(SP_2d{i})',100);
set(hh,'edgecolor','none');
%caxis([0 0.04]);
% caxis([-1 1]);
colormap(ax1,'default')
text(1.5,-2.5,'$d=10, t=188$','interpreter','latex',...
'fontsize',15,'color','w');
ylabel('Z','fontsize',12);
ylim([y{i}(1)+10-3.29 y{i}(1)+10+3.29]);
set(gca,'fontsize',12,'tickdir','out',...
    'ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;


ax2 = axes('position',[.35 .45 .26 .22]);
    [cc,hh] = contourf(x{i},y{i},squeeze(SP_sheardeform{i})',100);
set(hh,'edgecolor','none');
%caxis([-.002 .01]);
colormap(ax2,jet);
set(gca,'fontsize',12,'yticklabel',[],'tickdir','out',...
    'ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;
ylim([y{i}(1)+10-3.29 y{i}(1)+10+3.29]);


ax3 = axes('position',[.62 .45 .26 .22]);
  [cc,hh] = contourf(x{i},y{i},squeeze(BF3d{i})',100);set(hh,'edgecolor','none');
colormap(ax3,flipud((cbrewer('seq', 'YlGn', 100))));
caxis([-3.5,0]*1e-3);
set(gca,'fontsize',12,'tickdir','out','ticklength',[0.02 0.02],...
    'xminortick','on','yminortick','on');box on;
ylim([y{i}(1)+10-3.29 y{i}(1)+10+3.29]);
set(gca,'fontsize',12,'yticklabel',[],'tickdir','out',...
    'ticklength',[0.02 0.02],'xminortick','on','yminortick','on');

i=3;
 ax1 = axes('position',[.08 .13 .26 .22]);
[cc,hh] = contourf(x{i},y{i},squeeze(SP_2d{i}'),100);
set(hh,'edgecolor','none');
%caxis([0 0.03]);
colormap(ax1,'default')
ylabel('Z','fontsize',12);xlabel('X','fontsize',12);
col = colorbar('location','southoutside','position',[.09,0.056 0.19 0.015]);
set(gca,'fontsize',12,...
    'tickdir','out','ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;
text(1.5,-5.8,'$d=3, t=266$','interpreter','latex',...
'fontsize',15,'color','w');
ylim([y{i}(1) 0])
 
ax2 = axes('position',[.35 .13 .26 .22]);
[cc,hh] = contourf(x{i},y{i},squeeze(SP_sheardeform{i})',100);
set(hh,'edgecolor','none');
%caxis([-4e-3 0.01]);
colormap(ax2,jet);
xlabel('X','fontsize',12);
set(gca,'fontsize',12,'yticklabel',[],...
    'tickdir','out','ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;
ylim([y{i}(1) 0]);


ax3 = axes('position',[.62 .13 .26 .22]);
[cc,hh] = contourf(x{i},y{i},squeeze(BF3d{i})',100);set(hh,'edgecolor','none');
colormap(ax3,flipud((cbrewer('seq', 'YlGn', 100))));
col = colorbar('location','southoutside','position',[.75,0.056 0.19 0.015]);
xlabel('X','fontsize',12);
set(gca,'fontsize',12,'yticklabel',[],...
    'tickdir','out','ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
ylim([y{i}(1) 0]);








