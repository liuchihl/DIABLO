% This script shows how to load in 2D slices and make a movie of the simulation output
% Run after readmean.m, and is for boundary problems
addpath(genpath('/glade/u/home/liuchihl/TOOLS'));

load('/glade/u/home/liuchihl/work/colormaps/storm_color.mat');
% cm1 = load('/glade/u/home/liuchihl/work/colormaps/NCL_colormap/MPL_Set1.rgb');
% cm2 = load('/glade/u/home/liuchihl/work/colormaps/NCL_colormap/MPL_s3pcpn_l.rgb');
cm_eps = load('NCV_rainbow2.rgb');
cm_Bz = load('MPL_gnuplot.rgb');
d=2.5;
% base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_10/butterfly/Ri_0.16_0.05_9';
 base_dir=['/glade/campaign/univ/uosc0024/diablo_2.0/Large_case/d_6/boundary/Ri_0.12_0.1_3'];


 NY=361; LY=20;
 NX=512; LX=28.56;
%  NX=512; LX=39.27;
%  NX=256; LX=18.48;
% NX=256; LX=14.96;

% NY=541; LY=30;
% NX=512; LX=36.96;
i=1;
x=linspace(0,LX(i),NX(i)); y=gyf;%linspace(-LY(i)/2,LY(i)/2,NY(i)); 
% z=linspace(0,LZ(i),NZ(i));

%% animation
%close all;
clear plot_xy M
filename=[base_dir '/movie.h5'];
% plot_xy = VideoWriter('Boundary_0.12_6_3'.avi']);

%  plot_xy.FrameRate = 20; plot_xy.Quality = 100;
%  open(plot_xy);



%%

tt= 120:200;%:500%:280;%:500;
% tini=1
close all;

  figure('position',[500,500,500,500])
 for k=tt
%  figure('position',[50 50 1300 750])

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

varname1=['/th1_xy/' timename];
varname2=['/epsilon_xy/' timename];
varname3=['/u_xy/' timename];

TH=h5read(filename,varname1);
A2=mmderiv(gyf,TH')';
pcolor(x,gyf,TH');shading interp

caxis([-1 1]);
shading flat
z2 = colorbar('southoutside'); set(z2,'position',[.37,.06,0.3,.015]);
colormap(cm_eps/255);
%cmocean('ice')
text(2,12,'$B$','fontsize',18,'interpreter','latex');
xlabel('X','fontsize',12);%ylabel('Z','fontsize',12);
title(['Ri =' num2str(RI) ', d =' num2str(d) ', Time=' num2str(time(k),'%.1f')],'fontname','times','fontweight','normal')

text(2,12,'$\epsilon$','fontsize',25,'color','y','interpreter','latex');
set(gca,'tickdir','out','xminortick','on','yminortick','on','fontsize',12);
z3 = colorbar('southoutside'); set(z3,'position',[.68,.06,0.3,.015]);

clear M
%  if size(M.cdata,1)~=1306||size(M.cdata,2)~=680
% %      
% %  set(gcf,'units','centimeters','paperunits','centimeters')
% %  set(gcf,'PaperType','A4');
% %  pp=[0.63 1 38 19.8];
% %  ps=[0 0 pp(3)/1.1 pp(4)/1.1];
% %  set(gcf,'paperposition',pp)
% %  set(gcf,'position',ps)
%  drawnow
%   end
 drawnow
%  M=getframe(gcf);
%   writeVideo(plot_xy,M); 
 clf;
%  close
 end

%  close(plot_xy)

