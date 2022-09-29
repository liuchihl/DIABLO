
clear ; 
addpath(genpath('/glade/u/home/liuchihl/TOOLS'));
load('/glade/u/home/liuchihl/work/colormaps/storm_color.mat');
% cm1 = load('/glade/u/home/liuchihl/work/colormaps/NCL_colormap/MPL_Set1.rgb');
% cm2 = load('/glade/u/home/liuchihl/work/colormaps/NCL_colormap/MPL_s3pcpn_l.rgb');
cm_eps = load('NCV_rainbow2.rgb');
cm_Bz = load('MPL_gnuplot.rgb');
base_dir=['/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_3/boundary/Ri_0.12_0.1_9'];
% base_dir=['/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_10/butterfly/Ri_0.12_0.05_8'];
filename_mean=[base_dir '/mean.h5'];
    varname=['/gyf/' '0001'];             % Y-COORDINATE
    gyf=h5read(filename_mean,varname);

%     Ri=0.14
% NY=361; LY=20;
% NX=512; LX=28.08;
% NZ=128; LZ=7.02;
% Ri=0.12



% In DIABLO coordinate (Ri=0.16)
% NY=[721,541,397];         % Here, NY should match the value in grid_def.all
% LY = [19.86,14.75,13.16];     %vertical length
% NX=[1024,1024,1024];
% LX=[28.2 27.9 28.8];
% LZ=[7.05,6.98,7.2];  
% NZ=[256,256,256];

% cases with LY=20 
% NY=361; LY=20;
% NX=512; LX=28.28;
% NZ=128; LZ=7.07;


%  NY=361; LY=20;
%  NX=512; LX=28.21;
%  NZ=128; LZ=7.05;
% 
%d=3
RI=0.12; Re=1000; Pr=1;
LX=28.36; NX=512;
LY=20;    NY=361;
LZ=7.09;  NZ=128;

% 
% LX=29.16; NX=512;
% LY=20;    NY=361;
% LZ=7.29;  NZ=128;

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
% plot_xy = VideoWriter('Butterfly_12_8');

plot_xy = VideoWriter('Boundary_0.12_3_9.avi');
% 
plot_xy.FrameRate = 28; plot_xy.Quality = 100;
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
tt= 50:300%400;%100:2:380;%:1100%:500%:280;%:500;
% tt= 1052:1382;
% for k=1:501
%     if (k<10)
%         timename=['000' int2str(k)];
%     elseif (k<100)
%         timename=['00' int2str(k)];
%     elseif (k<1000)
%         timename=['0' int2str(k)];
%     else
%         timename=[int2str(k)];
%     end
% 
%     varname=['/time/' timename];            % TIME
%     time(k)=h5read(filename_mean,varname); 
% end
close all;
    figure;
set(gcf,'units','centimeters','paperunits','centimeters')
set(gcf,'PaperType','A4');
pp=[0.63 1 38 19.8];
ps=[0 0 pp(3)/1.1 pp(4)/1.1];
set(gcf,'paperposition',pp);
set(gcf,'position',ps);
 for k=tt
%  figure('position',[50 50 1300 750])
clf
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
varname2=['/epsilon_xy/' timename];
varname3=['/u_xy/' timename];
%  varname4=['/u_xz/' timename];
% varnamePlus = ['/uPlus_xz/' timename];
%  varnameMinus = ['/uMinus_xz/' timename];
 
TH=h5read(filename,varname1);
EPS=h5read(filename,varname2);
U=h5read(filename,varname3);
% U_xz=h5read(filename,varname4); 
% uPlus = h5read(filename,varnamePlus); 
% uMinus = h5read(filename,varnameMinus); 
% A1=mmderiv(gyf,U')';
% A2=mmderiv(gyf,TH')';
A3 = EPS; 
% [A3(k),po(k)] = max(mean(EPS));
%  
%   
ax1 = axes('position',[.06 .15 .3 .7]);

pcolor(x,gyf,U'); shading flat
caxis([-1 1]);
%  axis equal tight
colormap(ax1,jet);
% tit = title(sprintf('Time=%.2f',time(k)));
set(gca,'tickdir','out','xminortick','on','yminortick','on','fontsize',12);
z1 = colorbar('southoutside'); set(z1,'position',[.06,.06,0.3,.015]);
text(2,9,'$U$','fontsize',18,'interpreter','latex');
ylabel('$Z$','fontsize',12,'interpreter','latex');
% ylim([-5 5 ]);

ax2 = axes('position',[.37 .15 .3 .7]);
pcolor(x,gyf,TH');shading interp
%  caxis([-1 4]);
% caxis([0 2]);
caxis([-1 1]);
%  axis equal tight
shading flat
set(gca,'yticklabel',[],'tickdir','out','xminortick','on','yminortick','on','fontsize',12);
z2 = colorbar('southoutside'); set(z2,'position',[.37,.06,0.3,.015]);
colormap(ax2,flipud(cm))
%cmocean('ice')
text(2,9,'$B$','fontsize',18,'interpreter','latex');
xlabel('X','fontsize',12);%ylabel('Z','fontsize',12);
title(sprintf(['$Pr=1, Ri_0=0.12, d=3$, case 9, $t=%.2f$'],time(k)),'interpreter','latex')

ax3 = axes('position',[.68 .15 .3 .7]);
pcolor(x,gyf,log10(A3)'); shading flat
caxis([-6 -2]);
%  axis equal tight;
colormap(ax3,cm_eps/255);
text(2,9,'$\epsilon$','fontsize',25,'color','y','interpreter','latex');
% ylim([-5 5 ]);
set(gca,'yticklabel',[],'tickdir','out','xminortick','on','yminortick','on','fontsize',12);
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
M=getframe(gcf);
 writeVideo(plot_xy,M); 
% clf;
%  close
 end
% 
 close(plot_xy)

 

 
 %% phase speed
%  close all;
% [~,po] = min(abs(gyf+7.5));
%  figure; [cc,hh] = contourf(x,time(1:320),squeeze(U(:,po,1:320))',50); set(gca,'fontsize',15,'ticklength',[.02,.02]);
%  set(hh,'edgecolor','none');hold on;
%  [cc,hh] = contour(x,time(1:320),squeeze(U(:,po,1:320))',30,'k'); 
%  
%  xlabel('X'); ylabel('t'); grid minor
% title('d=2.5','fontsize',18)
%% saving data
tt = [1:501];
clear Bxz Byz
% tt=1:501;
base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_10/butterfly/Ri_0.14_0.05_8';
filename_mean=[base_dir '/mean.h5'];
    varname=['/gyf/' '0001'];             % Y-COORDINATE
%     gyf=h5read(filename_mean,varname);

    i=1;
x=linspace(0,LX(i),NX(i)); %y=gyf;%y=linspace(-LY(i)/2,LY(i)/2,NY(i)); 
z=linspace(0,LZ(i),NZ(i));
    Z=h5read(filename_mean,varname);
    X=x;
    Y=z;
filename=[base_dir '/movie.h5'];

 for k=1:length(tt)
%  figure('position',[50 50 1300 750])

 k
   if (tt(k)<10)
     timename=['000' int2str(tt(k))];
   elseif (tt(k)<100) 
     timename=['00' int2str(tt(k))];
   elseif (tt(k)<1000)
     timename=['0' int2str(tt(k))];
   else
     timename=[int2str(tt(k))];
   end
   
 varname1=['/th1_xy/' timename];
 Bxz(:,:,k)=h5read(filename,varname1);
%  varname2=['/th1_zy/' timename];
%  Byz(:,:,k)=h5read(filename,varname2);
%  varname3=['/th1_xz/' timename];
%  Bxy(:,:,k)=h5read(filename,varname3);
 
    varname=['/time/' timename];            % TIME
    time(k)=h5read(filename_mean,varname); 
 end
 t=time;
save B_xz_butterfly14_15 Bxz X Y Z t
%% plot the centered buoyancy evolution for poster (z=Lz/2 in DIABLO coordinate)
% base_dir='/glade/campaign/univ/uosc0024/diablo_2.0/Large_case/d_3/Ri_0.16';
% 
% filename=[base_dir '/movie.h5'];
% % tt=[400,600,697,800,1400];
% tt=[466,575,943,987,1300];
% % close all;
% 
% figure;
%  set(gcf,'units','centimeters','paperunits','centimeters')
%  set(gcf,'PaperType','A4');
%  pp=[0.63 1 27 8];
% ps=[0 0 pp(3)/1.1 pp(4)/1.1];
%  set(gcf,'paperposition',pp);
%  set(gcf,'position',ps);
% 
%  for k=tt
%  k
%    if (k<10)
%      timename=['000' int2str(k)];
%    elseif (k<100) 
%      timename=['00' int2str(k)];
%    elseif (k<1000)
%      timename=['0' int2str(k)];
%    else
%      timename=[int2str(k)];
%    end
%  var1=['/th1_xy/' timename];
%  TH1=h5read(filename,var1);
% %  var1=['/u_xy/' timename];
% %  var2=['/v_xy/' timename];
% %  var3=['/w_xy/' timename];
% %  U=h5read(filename,var1);
% %  V=h5read(filename,var2);
% %  W=h5read(filename,var3);
% 
%  if k==tt(1)
% ax1 = axes('position',[.06 .19 .18 .7]);
%  pcolor(x,y,TH1'); shading flat
%  caxis([0 2.0]);
% %  axis equal tight
% %  colormap(ax1,jet);
% %  cmocean('thermal');
%  % tit = title(sprintf('Time=%.2f',time(k)));
%  set(gca,'tickdir','out','xminortick','on','yminortick','on','fontsize',12);
% %  z1 = colorbar('southoutside'); set(z1,'position',[.06,.06,0.3,.015]);
% %  text(2,5,'$dU/dz$','fontsize',18,'interpreter','latex');
%  ylabel('$Z$','fontsize',12,'interpreter','latex');
% %  ylim([-4 4 ]);
% %  ylim([y(1)+5-4 y(1)+5+4 ]);
%  ylim([y(1) y(1)+3+5 ]);
% xlabel('$X$','fontsize',12,'interpreter','latex');
%  title(sprintf('$t=%.1f$',time(k)),'interpreter','latex')
% cmocean('thermal');
% 
%  elseif k==tt(2)
% ax2 = axes('position',[.245 .19 .18 .7]);
%  pcolor(x,y,TH1');shading interp
% %  caxis([-1 4]);
% caxis([0 2])
% %  axis equal tight
%  shading flat
%  set(gca,'yticklabel',[],'tickdir','out','xminortick','on','yminortick','on','fontsize',12);
%  colormap(ax2,flipud(cm))
%   title(sprintf('$t=%.1f$',time(k)),'interpreter','latex')
% %  ylim([-4 4 ]);
% %  ylim([y(1)+5-4 y(1)+5+4 ]);
%  ylim([y(1) y(1)+3+5 ]);
% cmocean('thermal');
% xlabel('$X$','fontsize',12,'interpreter','latex');
% 
%  elseif k==tt(3)
% ax3 = axes('position',[.435 .19 .18 .7]);
%  pcolor(x,y,TH1'); shading flat
% caxis([0 2])
% %  axis equal tight;
% %  colormap(ax3,cm_eps/255);
%  % ylim([-5 5 ]);
%  set(gca,'yticklabel',[],'tickdir','out','xminortick','on','yminortick','on','fontsize',12);
% %  z3 = colorbar('southoutside'); set(z3,'position',[.68,.06,0.3,.015]);
% %  ylim([-4 4 ]);
% %  ylim([y(1)+5-4 y(1)+5+4 ]);
%  ylim([y(1) y(1)+3+5 ]);
%   title(sprintf('$t=%.1f$',time(k)),'interpreter','latex')
% cmocean('thermal');
% xlabel('$X$','fontsize',12,'interpreter','latex');
% 
%  elseif k==tt(4)
% ax4 = axes('position',[.625 .19 .18 .7]);
%  pcolor(x,y,TH1'); shading flat
% caxis([0 2])
% %  axis equal tight;
% %  colormap(cm_eps/255);
%  set(gca,'yticklabel',[],'tickdir','out','xminortick','on','yminortick','on','fontsize',12);
% %  z3 = colorbar('southoutside'); set(z3,'position',[.68,.06,0.3,.015]);
% %  ylim([-4 4 ]);
% %  ylim([y(1)+5-4 y(1)+5+4 ]);
%  ylim([y(1) y(1)+3+5 ]);
%   title(sprintf('$t=%.1f$',time(k)),'interpreter','latex')
% cmocean('thermal');
% xlabel('$X$','fontsize',12,'interpreter','latex');
% 
%  else 
% ax5 = axes('position',[.815 .19 .18 .7]);
%  pcolor(x,y,TH1'); shading flat
% caxis([0 2])
% %  axis equal tight;
%  set(gca,'yticklabel',[],'tickdir','out','xminortick','on','yminortick','on','fontsize',12);
%  z5 = colorbar('southoutside'); set(z5,'position',[.815,.08,0.18,.015]);
% %  ylim([-4 4 ]);
% %  ylim([y(1)+5-4 y(1)+5+4 ]);
%  ylim([y(1) y(1)+3+5 ]);
% 
%   title(sprintf('$t=%.1f$',time(k)),'interpreter','latex')
% cmocean('thermal');
% % xlabel('$X$','fontsize',12,'interpreter','latex');
% 
%  end
% %  colormap(cm_eps/255);
% 
%  end
 
%  print -djpeg -r300 B_evolution_3.jpg

%% for poster: pairing, draining mode B and epsilon
% i=1;
% % close all;
% filename=[base_dir '/movie.h5'];
% 
%      figure;
%  set(gcf,'units','centimeters','paperunits','centimeters')
%  set(gcf,'PaperType','A4');
%  pp=[0.63 1 27 10];
% ps=[0 0 pp(3)/1.1 pp(4)/1.1];
%  set(gcf,'paperposition',pp);
%  set(gcf,'position',ps);
%  tt=[929,923,1000];
%  k=tt(i)
% %  figure('position',[50 50 1300 750])
% %  ax1 = axes('position',[.06 .15 .4 .7]);
% %  ax2 = axes('position',[.47 .15 .4 .7]);
% 
% 
%    if (k<10)
%      timename=['000' int2str(k)];
%    elseif (k<100) 
%      timename=['00' int2str(k)];
%    elseif (k<1000)
%      timename=['0' int2str(k)];
%    else
%      timename=[int2str(k)];
%    end
%  
%  varname1=['/th1_xy/' timename];
%  varname2=['/epsilon_xy/' timename];
% %  var1=['/th1_xy/' timename];
% %  TH1=h5read(filename,var1);
% 
%  TH=h5read(filename,varname1);
%  EPS=h5read(filename,varname2);
% %  A1=mmderiv(gyf,U')';
% %  A2=mmderiv(gyf,TH')';
% %  A3=EPS;
%  
% ax1 = axes('position',[.065 .26 .45 .7]);
% % 
%  pcolor(x,y,TH'); shading flat
%  caxis([0 2]);
% %  axis equal tight
%  cmocean('thermal');
%  % tit = title(sprintf('Time=%.2f',time(k)));
%  set(gca,'tickdir','out','xminortick','on','yminortick','on','fontsize',12);
% %  z1 = colorbar('southoutside'); set(z1,'position',[.065,.1,0.45,.015]);
%  text(2,-5,'$B$','fontsize',18,'interpreter','latex');
%  ylabel('$Z$','fontsize',12,'interpreter','latex');
% %   xlabel('$X$','fontsize',12,'interpreter','latex');
% 
% %  ylim([-6 6 ]);
%  ylim([y(1) y(1)+12]);
%  ax2 = axes('position',[.525 .26 .45 .7]);
%  pcolor(x,y,log10(EPS)'); shading flat
%  caxis([-6 -2]);
% %  axis equal tight;
%  colormap(ax2,cm_eps/255);
%  text(2,-5,'$\epsilon$','fontsize',25,'color','y','interpreter','latex');
% %  ylim([-6 6 ]);
%  ylim([y(1) y(1)+12]);
%  set(gca,'yticklabel',[],'tickdir','out','xminortick','on','yminortick','on','fontsize',12);
% %  z2 = colorbar('southoutside'); set(z2,'position',[.525,.1,0.45,.015]);
% %   xlabel('$X$','fontsize',12,'interpreter','latex');
% 
% 
% 
% % print -djpeg -r300 B_eps_5.jpg




