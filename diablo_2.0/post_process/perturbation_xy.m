% This script shows how to load in 2D slices and make a movie of the simulation output
% Run after readmean.m

addpath(genpath('/glade/u/home/liuchihl/TOOLS'));
load('/glade/u/home/liuchihl/work/colormaps/storm_color.mat');
% cm1 = load('/glade/u/home/liuchihl/work/colormaps/NCL_colormap/MPL_Set1.rgb');
% cm2 = load('/glade/u/home/liuchihl/work/colormaps/NCL_colormap/MPL_s3pcpn_l.rgb');
cm_eps = load('NCV_rainbow2.rgb');
cm_Bz = load('MPL_gnuplot.rgb');
for cases = [12]
base_dir=['/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_10/butterfly/Ri_0.14_0.05_',num2str(cases)];
filename_mean=[base_dir '/mean.h5'];
    varname=['/gyf/' '0001'];             % Y-COORDINATE
    gyf=h5read(filename_mean,varname);

% Ri=0.12

% NY=[1009,1081,1081];         % Here, NY should match the value in grid_def.all
% LY = [28.2,30.375,31.4];     %vertical length
% NX=[1008,1024,1024];
% LX=[28.2 28.8 31.4];
% LZ=[7,7.2,7.9];  
% NZ=[252,256,245];


% In DIABLO coordinate (Ri=0.16)
% NY=[721,541,397];         % Here, NY should match the value in grid_def.all
% LY = [19.86,14.75,13.16];     %vertical length
% NX=[1024,1024,1024];
% LX=[28.2 27.9 28.8];
% LZ=[7.05,6.98,7.2];  
% NZ=[256,256,256];

% cases with LY=20 
% NY=721; LY=20;
% NX=1024; LX=28.61;
% NZ=256; LZ=7.15;


 NY=361; LY=20;
 NX=512; LX=28.08;
 NZ=128; LZ=7.02;

% test (small)
%fname1= 'KH_test2';
%LX=30; LY=20; LZ=4;
%NX=128;NY=65; NZ=16;
% timestep=[100:100:4000];
% after determining the grid sizes and domain size, compute x,y,z

i=1;
x=linspace(0,LX(i),NX(i)); y=gyf;%y=linspace(-LY(i)/2,LY(i)/2,NY(i)); 
z=linspace(0,LZ(i),NZ(i));

%% animation
%close all;
clear plot_xy M
filename=[base_dir '/movie.h5'];
% plot_xy = VideoWriter('Butterfly9');

plot_xy = VideoWriter(['Butterfly_14_perturbation_' num2str(cases) '.avi']);

plot_xy.FrameRate = 20; plot_xy.Quality = 100;
open(plot_xy);

% tt14=time([1:200,211:411,413:604,623:705,707:905,907:1079,1081:nk]);
% tt5=[1:159,169:360,362:611,613:862,864:1113,1115:nk];
% tt=[310:604,623:nk];
%            figure('position',[20,20,1200 700]); 
% figure; plot(diff((tt5)));
% tt5=[330:360,362:611,613:862,864:1113,1115:nk];
% tt10=[600:702,704:925,949:1176,1178:1552,1554:nk(1)];
% tt5=[600:852,854:1303,1305:nk];
%%
tt= 1:320;%100:2:380;%:1100%:500%:280;%:500;
% tini=1
close all;
% tt3=[400:753,755:1004,1006:1255,1257:nk];
% close all;figure;
% tt=[tt5];
%  tt3=nk;
% clear plot_xy
% ylim([-10, 3])

    figure;
set(gcf,'units','centimeters','paperunits','centimeters')
set(gcf,'PaperType','A4');
pp=[0.63 1 38 19.8];
ps=[0 0 pp(3)/1.1 pp(4)/1.1];
set(gcf,'paperposition',pp);
set(gcf,'position',ps);
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
    varname=['/time/' timename];            % TIME
    time(k)=h5read(filename_mean,varname); 
% varname = ['/ume_z/' timename];
% UME_Z = h5read(filename,varname); 
varname1=['/th1_xy/' timename];
varname3=['/u_xy/' timename];
 
TH=h5read(filename,varname1);
U=h5read(filename,varname3);
th = TH-mean(TH,1);
u = U-mean(U,1);
uz=mmderiv(gyf,u')';
 
  
ax1 = axes('position',[.06 .15 .4 .7]);

pcolor(x,gyf,uz'); shading flat
% caxis([-1 1]);
%  axis equal tight
colormap(ax1,jet);
% tit = title(sprintf('Time=%.2f',time(k)));
set(gca,'tickdir','out','xminortick','on','yminortick','on','fontsize',12);
z1 = colorbar('southoutside'); set(z1,'position',[.06,.06,0.4,.015]);
text(2,5,'$du''/dz$','fontsize',18,'interpreter','latex');
ylabel('$Z$','fontsize',12,'interpreter','latex');
% ylim([-5 5 ]);

ax2 = axes('position',[.49 .15 .4 .7]);
pcolor(x,gyf,th');shading interp
%  caxis([-1 4]);
% caxis([0 2]);
% caxis([-1 1]);
%  axis equal tight
shading flat
set(gca,'yticklabel',[],'tickdir','out','xminortick','on','yminortick','on','fontsize',12);
z2 = colorbar('southoutside'); set(z2,'position',[.49,.06,0.4,.015]);
colormap(ax2,flipud(cm))
%cmocean('ice')
text(2,5,'$b''$','fontsize',18,'interpreter','latex');
xlabel('X','fontsize',12);%ylabel('Z','fontsize',12);
title(sprintf(['$Pr=1, Ri_0=0.14, case ' num2str(cases) ' , Time=%.2f$'],time(k)),'interpreter','latex')


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
M=getframe(gcf);
 writeVideo(plot_xy,M); 
clf;
%  close
 end
% 
 close(plot_xy)


end
