% Run after readmean.m
addpath(genpath('/glade/u/home/liuchihl/TOOLS'));

load('/glade/u/home/liuchihl/work/colormaps/storm_color.mat');
% cm1 = load('/glade/u/home/liuchihl/work/colormaps/NCL_colormap/MPL_Set1.rgb');
% cm2 = load('/glade/u/home/liuchihl/work/colormaps/NCL_colormap/MPL_s3pcpn_l.rgb');
cm_eps = load('NCV_rainbow2.rgb');
cm_Bz = load('MPL_gnuplot.rgb');
 
LZ=7.02;
NZ=128;
% load('/glade/u/home/liuchihl/work/colormaps/storm_color.mat')
cm1 = load('/glade/u/home/liuchihl/work/colormaps/NCL_colormap/MPL_Set1.rgb');
cm2 = load('/glade/u/home/liuchihl/work/colormaps/NCL_colormap/MPL_s3pcpn_l.rgb');
cm_Bz = load('MPL_gnuplot.rgb');

% Background density gradient
drhodz1=0.0;
z=linspace(0,LZ,NZ);    % spanwise direction
filename=[base_dir '/movie.h5'];
%%
clear plot_yz 
plot_yz = VideoWriter('Butterfly_14_12yz.avi');
plot_yz.FrameRate = 20; plot_yz.Quality = 100;
open(plot_yz);
tt3 = 50:501;
close all;

 figure;
 set(gcf,'units','centimeters','paperunits','centimeters')
 set(gcf,'PaperType','A4');
 pp=[0.63 1 38 19.8];
ps=[0 0 pp(3)/1.1 pp(4)/1.1];
 set(gcf,'paperposition',pp);
 set(gcf,'position',ps);
for k=tt3

k

  if (k<10)
    timename=['000' int2str(k)];
  elseif (k<100) 
    timename=['00' int2str(k)];
  elseif (k<1000)
    timename=['0' int2str(k)];
  else
    timename=[int2str(k)];
  end

 varname1=['/th1_zy/' timename];
 varname2=['/epsilon_zy/' timename];
  varname4=['/u_zy/' timename];

 TH=h5read(filename,varname1);
 EPS=h5read(filename,varname2);
 U=h5read(filename,varname4); 
 A1=mmderiv(gyf,U')';
 A2=mmderiv(gyf,TH')';
 A3=EPS;
 
% surf(z,gyf,zeros(size(A_nu_t')),A_nu_t','EdgeColor','none'),view(0,90);
% hold on
% caxis([0 1e-4]);
% contour(z/1e3,gyf-gyf(end),A_v',10);
ax1 = axes('position',[.06 .15 .3 .7]);
pcolor(z,gyf,A1');shading interp;
caxis([-1 1.5]); colormap(ax1,jet);
set(gca,'tickdir','out','xminortick','on','yminortick','on','fontsize',12);
z1 = colorbar('southoutside'); set(z1,'position',[.13,.06,0.15,.015]);
% ylim([-4,-1]);
text(1.1,-6,'$dU/dz$','fontsize',18,'interpreter','latex');
ylabel('$Z$','fontsize',12,'interpreter','latex');
xlabel('$Y$','fontsize',12,'interpreter','latex');%ylabel('Z','fontsize',12);

ax2 = axes('position',[.37 .15 .3 .7]);
pcolor(z,gyf,TH');shading interp
 caxis([-1 1]);
% axis equal tight
% ylim([gyf(1) gyf(1)+5])
% ylim([-4,-1]);
 set(gca,'yticklabel',[],'tickdir','out','xminortick','on','yminortick','on','fontsize',12);
 z2 = colorbar('southoutside'); set(z2,'position',[.44,.06,0.15,.015]);
 colormap(ax2,flipud(cm))
 %cmocean('ice')
%  text(1.1,-6,'$B$','fontsize',18,'interpreter','latex');
 xlabel('$Y$','fontsize',12,'interpreter','latex');%ylabel('Z','fontsize',12);
 title(sprintf('$Pr=1, Ri_0=0.14, case12, Time=%.2f$',time(k)),'interpreter','latex')
 
 ax3 = axes('position',[.68 .15 .3 .7]);
 pcolor(z,gyf,log10(A3)'); shading interp
 caxis([-6 -3]);
%  axis equal tight;
%  ylim([gyf(1) gyf(1)+5])
% ylim([-4,-1]);
 colormap(ax3,cm_eps/255);
 text(1.1,-6,'$\epsilon$','fontsize',25,'color','y','interpreter','latex');
 % ylim([-5 5 ]);
 xlabel('$Y$','fontsize',12,'interpreter','latex');%ylabel('Z','fontsize',12);
 set(gca,'yticklabel',[],'tickdir','out','xminortick','on','yminortick','on','fontsize',12);
 z3 = colorbar('southoutside'); set(z3,'position',[.75,.06,0.15,.015]);
clear M;

 set(gcf,'units','centimeters','paperunits','centimeters')
 set(gcf,'PaperType','A4');
 pp=[0.63 1 38 19.8];
ps=[0 0 pp(3)/1.1 pp(4)/1.1];
 set(gcf,'paperposition',pp);
 set(gcf,'position',ps);

M=getframe(gcf);
clf;
% close
writeVideo(plot_yz,M);
end
close(plot_yz)
%print -djpeg coarse_d=3_Lz4.jpg
