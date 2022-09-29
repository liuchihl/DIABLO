% This script illustrates how to load in 3D model data and calculate 
% some basic diagnostics
% Run after readmean.m
clear;%clc; 

% input coefficients
fname1 = 'd_10/Ri_0.16/';  
fname2= 'd_5/Ri_0.16/';  
fname3 = 'd_3/Ri_0.16/';  
% In DIABLO coordinate
NY=[721,541,397];         % Here, NY should match the value in grid_def.all
LY = [19.86,14.75,13.16];     %vertical length
NX=1024;
LX=[28.2 27.9 28.8];
LZ=[7.05,6.98,7.2];  
NZ=256;
RI=0.16; Re=1000; Pr=1;
% timestep = [22000,29000;23000,30000;32000,40000];
% timestep = [22000;23000;32000];
timestep = [26000;27000;36000];
% timestep = [28000;29000;38000];

%  timestep = [1000:1000:60000];
% peak of APE
% d=10: t=153 (22000), 200 (29000)
% d=5: t=154 (23000), 201 (30000)
% d=3: t=221 (32000), 277 (40000)
% timestep = [36000];
addpath(genpath('/glade/u/home/liuchihl/TOOLS'));
cm_1 = load('MPL_gnuplot.rgb');cm_eps = load('NCV_rainbow2.rgb');

for i = 1:3
U=zeros(NX,NY(i),NZ);V=zeros(NX,NY(i),NZ);
W=zeros(NX,NY(i),NZ);TH1=zeros(NX,NY(i),NZ);
U3d=zeros(NX,NY(i),NZ,length(timestep));
V3d=zeros(NX,NY(i),NZ,length(timestep));
W3d=zeros(NX,NY(i),NZ,length(timestep));TH3d=zeros(NX,NY(i),NZ,length(timestep));
U2d=zeros(NX,NY(i),NZ);
% V2d=zeros(NX,NY(i),NZ);W2d=zeros(NX,NY(i),NZ);
% U3d_x=zeros(NX,NY(i),NZ);V3d_y=zeros(NX,NY(i),NZ);
% V3d_x=zeros(NX,NY(i),NZ);U3d_y=zeros(NX,NY(i),NZ);
%W3d_z=zeros(NX,NY(i),NZ);%SP_background{i}=zeros(size(timestep));
% SP_2d{i}=zeros(size(timestep));SP_background{i}=zeros(size(timestep));
% SD3d{i}=zeros(size(timestep));BF3d{i}=zeros(size(timestep));
% D3d{i}=zeros(size(timestep)); SP_sheardeform{i}=zeros(size(timestep));
TH1y{i}=zeros(NX,NY(i),NZ);
TH3d_sq{i}=zeros(NX,NY(i));
K3D{i}=zeros(NX,NY(i),length(timestep(i,:)));
TH3d_x=zeros(NX,NY(i),NZ); TH3d_y=zeros(NX,NY(i),NZ);
TH3d_z=zeros(NX,NY(i),NZ);
BF3d_v{i}=zeros(NX,NY(i)); BF3d_h{i}=zeros(NX,NY(i));
dTHzdy{i}=zeros(NX,NY(i)); dTHzdx{i}=zeros(NX,NY(i));
eps3d_p{i}=zeros(size(timestep));

if     i==1    
    base_dir=['/glade/campaign/univ/uosc0024/diablo_2.0/Large_case/',fname1,'output/'];
elseif i==2
    base_dir=['/glade/campaign/univ/uosc0024/diablo_2.0/Large_case/',fname2,'output/'];
else
    base_dir=['/glade/campaign/univ/uosc0024/diablo_2.0/Large_case/',fname3,'output/'];

end
filename=[base_dir];
for k=1:length(timestep(i,:))

  if (timestep(i,k)<10)
    timename=['out0000' int2str(timestep(i,k)) '.h5'];
  elseif (timestep(i,k)<100)
    timename=['out000' int2str(timestep(i,k)) '.h5'];
  elseif (timestep(i,k)<1000)
    timename=['out00' int2str(timestep(i,k)) '.h5'];
  elseif (timestep(i,k)<10000)
    timename=['out0' int2str(timestep(i,k)) '.h5'];
  else
    timename=['out' int2str(timestep(i,k)) '.h5'];
  end
varname1=['/Timestep/U'];
varname2=['/Timestep/V'];
varname3=['/Timestep/W'];
varname4=['/Timestep/TH1'];
% index for time
% if k==1
%    time = 1;
%  elseif k==timestep(2)
%      time = 2;
% % else k==timestep(3)
% %     time = 3;
% end
if exist([filename,timename])
U(:,:,:)=h5read([filename,timename],varname1);
V(:,:,:)=h5read([filename,timename],varname2);
W(:,:,:)=h5read([filename,timename],varname3);
TH1(:,:,:)=h5read([filename,timename],varname4);

info = h5info([filename,timename]);
%time
tt{i}(k) = info.Groups.Attributes.Value;
% NX=size(U,1); NZ=size(U,3); NY=size(U,2);
x{i}=linspace(0,LX(i),NX); y{i}=linspace(-LY(i)/2,LY(i)/2,NY(i)); 
z{i}=linspace(0,LZ(i),NZ);

zz = NZ/2;

% calculate 3D perturbation kinetic energy (Smyth et al., 2005)
U3d = U-mean(U,3);V3d = V-mean(V,3);W3d = W-mean(W,3);
TH3d = TH1-mean(TH1,3);
K3d{i}(k) = mean(mean(mean(.5*(U3d.^2+V3d.^2+W3d.^2))));  %

for j=1:NY(i)
    TH3d_x(:,j,:)= gradient(squeeze(TH3d(:,j,:))',mean(diff(x{i})))';
end
for j=1:NX
    TH3d_y(j,:,:) = gradient(squeeze(TH3d(j,:,:))',mean(diff(y{i})))';
    TH3d_z(j,:,:) = gradient(squeeze(TH3d(j,:,:)),mean(diff(z{i})));
end
% U_mean = mean(U,3); U_rms = sqrt((U-U_mean).^2);
% V_mean = mean(V,3); V_rms = sqrt((V-V_mean).^2);
% W_mean = mean(W,3); W_rms = sqrt((W-W_mean).^2);
end


TH3d_sq{i}(:,:) = mean(TH3d.^2,3);        %spanwise-averaged TH3d^2

% calculate the buoyancy variance budget equation
% dissipation rate (epsilon_p) 
dTHzdy{i} = gradient(squeeze(mean(TH1,3)),mean(diff(y{i})));  
dTHzdx{i} = gradient(squeeze(mean(TH1,3))',mean(diff(x{i})))';   

eps3d_p{i} = -1/(Re*Pr)*mean(...
    (TH3d_x.^2+TH3d_y.^2+TH3d_z.^2),3);
% buoyancy flux
BF3d_v{i} = -mean(dTHzdy{i}.*TH3d.*V3d,3);
BF3d_h{i} = -mean(dTHzdx{i}.*TH3d.*U3d,3);



% K3D{i}(:,:,k) = squeeze(mean(.5*(U3d(:,:,:).^2+V3d(:,:,:).^2+...
%     W3d(:,:,:).^2),3));



end

end
%tke_3D_m{i}
%SP_background{i}+SP_2d{i}+SP_sheardeform{i}+BF3d{i}+D3d{i}
%  save K3D_budget_d3 SP_background SP_2d SD3d SP_sheardeform BF3d D3d tke_3D_m tke_2D_m tt
% save('omega_K3D_16_3_20000_30000','omega_x','K3D','tt','-v7.3')


% plot the averaged omega_x^2 and K3D in x-z plane
% omega_x_sq_mean = mean(omega_x.^2,3);
% omega_y_sq_mean = mean(omega_y.^2,3);
% omega_z_sq_mean = mean(omega_z.^2,3);
% save center_omega_K3d K3D TH1y_cen omega_x_cen omega_y_cen omega_z_cen tt
% save center_omega_K3d K3D TH1y_cen omega_x_cen omega_y_cen omega_z_cen tt
% save center_omega_K3d_onset K3D TH1y_cen omega_x_cen omega_y_cen omega_z_cen tt





%% %% plot the x-z of the buoyancy variance budget  TH3d^2/2

close all;
figure;
 set(gcf,'units','centimeters','paperunits','centimeters')
 set(gcf,'PaperType','A4');
 pp=[0.63 .9 26 24];
ps=[0 0 pp(3)/1.1 pp(4)/1.1];
 set(gcf,'paperposition',pp);
 set(gcf,'position',ps);

 
 
%  k=1;%:length(timestep)
i=1;

ax1 = axes('position',[.08 .725 .21 .22]);
  [cc,hh] = contourf(x{i},y{i},squeeze(.5*TH3d_sq{i})',100);
set(hh,'edgecolor','none');
title(' $b^2_{3d}/2$','fontname','times','fontsize',18,...
    'fontweight','normal','interpreter','latex');
caxis([0 0.04]);
% caxis([-1 1]);
colormap(ax1,'default')

% hold on; contour(x{i},y{i},TH1y_cen{i}(:,:,k)',[.6,1],'k');
% colormap(ax1,(cbrewer('seq', 'PuBuGn', 100)));
% colorbar('position',[0.85,0.55 0.015 0.4]);
text(1.5,-2.5,'$d=10, t=188$','interpreter','latex',...
'fontsize',15,'color','w');
ylabel('Z','fontsize',12);
ylim([y{i}(1)+10-3.29 y{i}(1)+10+3.29]);
set(gca,'fontsize',12,'tickdir','out',...
    'ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;
% title('$Re$=1000, $Pr$=1, $Ri$=0.16, $d$=3','fontname','times','fontsize',12,...
%     'fontweight','normal','interpreter','latex');
% ylim([y{i}(1) 0]);

  ax2 = axes('position',[.3 .725 .21 .22]);
  [cc,hh] = contourf(x{i},y{i},squeeze(BF3d_h{i})',100);
set(hh,'edgecolor','none');
caxis([-.002 .01]);
% hold on; contour(x{i},y{i},TH1y_cen{i}(:,:,k)',[.6,1],'k');
% colormap(ax2,flipud(cbrewer('seq', 'YlOrRd', 100)));
colormap(ax2,jet);
% colorbar('position',[0.85,0.55 0.015 0.4]);
title('$BF_{3d}^{hor}$','interpreter','latex',...
'fontsize',18,'color','k');
% xlabel('X'); 
% ylabel('Z');
set(gca,'fontsize',12,'yticklabel',[],'tickdir','out',...
    'ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;
ylim([y{i}(1)+10-3.29 y{i}(1)+10+3.29]);


ax3 = axes('position',[.52 .725 .21 .22]);
  [cc,hh] = contourf(x{i},y{i},squeeze(BF3d_v{i})',100);
set(hh,'edgecolor','none');
caxis([-.002 .01]);
% hold on; contour(x{i},y{i},TH1y_cen{i}(:,:,k)',[.6,1],'k');
colormap(ax3,jet);

% colormap(ax3,flipud(cbrewer('seq', 'YlOrRd', 100)));
% colorbar('position',[0.85,0.55 0.015 0.4]);
title('$BF_{3d}^{ver}$','interpreter','latex',...
'fontsize',18,'color','k');
% xlabel('X'); 
% ylabel('Z');
set(gca,'fontsize',12,'yticklabel',[],'tickdir','out',...
    'ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;
% ylim([y{i}(1) 0])
ylim([y{i}(1)+10-3.29 y{i}(1)+10+3.29]);


  ax4 = axes('position',[.74 .725 .21 .22]);
  [cc,hh] = contourf(x{i},y{i},squeeze(eps3d_p{i})',100);set(hh,'edgecolor','none');
% hold on; [cc,hh]=contour(x{i},y{i},tke_3D_m',10,'color','k');%set(cc,'color','k');
colormap(ax4,flipud((cbrewer('seq', 'YlGn', 100))));
caxis([-3.5,0]*1e-3);
% hold on; contour(x{i},y{i},TH1y_cen{i}(:,:,k)',[.6 1],'k');
% colorbar('position',[0.85,0.12 0.015 0.4]);
title('$\epsilon_{p_{3d}}$','interpreter','latex',...
'fontsize',15,'color','k');
% xlabel('X'); ylabel('Z');
set(gca,'fontsize',12,'tickdir','out','ticklength',[0.02 0.02],...
    'xminortick','on','yminortick','on');box on;
% ylim([y{i}(1) 0]);
ylim([y{i}(1)+10-3.29 y{i}(1)+10+3.29]);
set(gca,'fontsize',12,'yticklabel',[],'tickdir','out',...
    'ticklength',[0.02 0.02],'xminortick','on','yminortick','on');

i=2;
  ax1 = axes('position',[.08 .45 .21 .22]);
  [cc,hh] = contourf(x{i},y{i},squeeze(.5*TH3d_sq{i})',100);
set(hh,'edgecolor','none');
caxis([0 0.04]);
col = colorbar('location','southoutside','position',[.09,0.4 0.19 0.015]);
colormap(ax1,'default')

% hold on; contour(x{i},y{i},TH1y_cen{i}(:,:,k)',[.6,1],'k');
% colormap(ax1,(cbrewer('seq', 'PuBuGn', 100)));
% colorbar('position',[0.85,0.55 0.015 0.4]);
text(1.5,-5,'$d=5, t=193$','interpreter','latex',...
'fontsize',15,'color','w');
ylabel('Z','fontsize',12);
ylim([y{i}(1)+5-3.29 y{i}(1)+5+3.29]);
set(gca,'fontsize',12,'tickdir','out',...
    'ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;
% title('$Re$=1000, $Pr$=1, $Ri$=0.16, $d$=3','fontname','times','fontsize',12,...
%     'fontweight','normal','interpreter','latex');
% ylim([y{i}(1) 0]);

  ax2 = axes('position',[.3 .45 .21 .22]);
  [cc,hh] = contourf(x{i},y{i},squeeze(BF3d_h{i})',100);
set(hh,'edgecolor','none');
caxis([-.002 .01]);
colormap(ax2,jet);
% hold on; contour(x{i},y{i},TH1y_cen{i}(:,:,k)',[.6,1],'k');
% colormap(ax2,flipud(cbrewer('seq', 'YlOrRd', 100)));
% col = colorbar('location','southoutside','position',[.31,0.37 0.19 0.015]);
% xlabel('X'); 
% ylabel('Z');
set(gca,'fontsize',12,'yticklabel',[],'tickdir','out',...
    'ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;
ylim([y{i}(1)+5-3.29 y{i}(1)+5+3.29]);


ax3 = axes('position',[.52 .45 .21 .22]);
  [cc,hh] = contourf(x{i},y{i},squeeze(BF3d_v{i})',100);
set(hh,'edgecolor','none');
caxis([-.002 .01]);
% hold on; contour(x{i},y{i},TH1y_cen{i}(:,:,k)',[.6,1],'k');
colormap(ax3,jet);

% colormap(ax3,flipud(cbrewer('seq', 'YlOrRd', 100)));
col = colorbar('location','southoutside','position',[.53,0.4 0.19 0.015]);
set(col,'xtick',[0:0.01:0.02],'xticklabel',{'0','0.01','0.02'});
% xlabel('X'); 
% ylabel('Z');
set(gca,'fontsize',12,'yticklabel',[],'tickdir','out',...
    'ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;
% ylim([y{i}(1) 0])
ylim([y{i}(1)+5-3.29 y{i}(1)+5+3.29]);

ax4 = axes('position',[.74 .45 .21 .22]);
 [cc,hh] = contourf(x{i},y{i},squeeze(eps3d_p{i})',100);set(hh,'edgecolor','none');
% hold on; [cc,hh]=contour(x{i},y{i},tke_3D_m',10,'color','k');%set(cc,'color','k');
colormap(ax4,flipud((cbrewer('seq', 'YlGn', 100))));
caxis([-3.5,0]*1e-3);
col = colorbar('location','southoutside','position',[.75,0.4 0.19 0.015]);

set(gca,'fontsize',12,'tickdir','out','ticklength',[0.02 0.02],...
    'xminortick','on','yminortick','on');box on;
% ylim([y{i}(1) 0]);
ylim([y{i}(1)+5-3.29 y{i}(1)+5+3.29]);
set(gca,'fontsize',12,'yticklabel',[],'tickdir','out',...
    'ticklength',[0.02 0.02],'xminortick','on','yminortick','on');


i=3;

  ax1 = axes('position',[.08 .13 .21 .22]);
[cc,hh] = contourf(x{i},y{i},squeeze(.5*TH3d_sq{i})',100);
set(hh,'edgecolor','none');
caxis([0 0.03]);
colormap(ax1,'default')
% hold on; contour(x{i},y{i},TH1y_cen{i}(:,:,k)',[.6,1],'k');
% colormap(ax1,(cbrewer('seq', 'PuBuGn', 100)));
% text(1.5,y{i}(1)*.8,'$\langle \omega_{x}^{2} \rangle_{y}$','interpreter','latex',...
% 'fontsize',12,'color','w');
ylabel('Z','fontsize',12);xlabel('X','fontsize',12);
col = colorbar('location','southoutside','position',[.09,0.056 0.19 0.015]);

set(gca,'fontsize',12,...
    'tickdir','out','ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;
text(1.5,-5.8,'$d=3, t=266$','interpreter','latex',...
'fontsize',15,'color','w');
ylim([y{i}(1) 0])

  ax2 = axes('position',[.3 .13 .21 .22]);

  [cc,hh] = contourf(x{i},y{i},squeeze(BF3d_h{i})',100);
set(hh,'edgecolor','none');
caxis([-4e-3 0.01]);
colormap(ax2,jet);
% hold on; contour(x{i},y{i},TH1y_cen{i}(:,:,k)',[.6,1],'k');
% colorbar('position',[0.9,0.53 0.015 0.2]);
% text(1.5,y{i}(1)*.8,'$\langle \omega_{y}^{2} \rangle_{y}$','interpreter','latex',...
% 'fontsize',12,'color','w');
xlabel('X','fontsize',12);
% ylabel('Z');
set(gca,'fontsize',12,'yticklabel',[],...
    'tickdir','out','ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;
% title('$Re$=1000, $Pr$=1, $Ri$=0.16, $d$=3','fontname','times','fontsize',12,...
%     'fontweight','normal','interpreter','latex');
ylim([y{i}(1) 0])

  ax3 = axes('position',[.52 .13 .21 .22]);
[cc,hh] = contourf(x{i},y{i},squeeze(BF3d_v{i})',100);
set(hh,'edgecolor','none');
caxis([-4e-3 0.01]);
colormap(ax3,jet);
col = colorbar('location','southoutside','position',[.53,0.056 0.19 0.015]);
% set(col,'ytick',[-.02:.01:.02])
% text(1.5,y{i}(1)*.8,'$\langle \omega_{z}^{2} \rangle_{y}$','interpreter','latex',...
% 'fontsize',12,'color','w');
% xlabel('X'); 
% ylabel('Z');
set(gca,'fontsize',12,'yticklabel',[],...
    'tickdir','out','ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;
xlabel('X','fontsize',12);
ylim([y{i}(1) 0])

ax4 = axes('position',[.74 .13 .21 .22]);
[cc,hh] = contourf(x{i},y{i},squeeze(eps3d_p{i})',100);set(hh,'edgecolor','none');
% hold on; [cc,hh]=contour(x{i},y{i},tke_3D_m',10,'color','k');%set(cc,'color','k');
colormap(ax4,flipud((cbrewer('seq', 'YlGn', 100))));
caxis([-3.5,0]*1e-3);
% hold on; contour(x{i},y{i},TH1y_cen{i}(:,:,k)',[.6,1],'k');
col = colorbar('location','southoutside','position',[.75,0.056 0.19 0.015]);
% set(col,'ytick',[0:5e-4:1.5e-3],'yticklabel',{'0','0.0005','0.001','0.0015'})

% text(1.5,y{i}(1)*.8,'$K_{3D}$','interpreter','latex',...
% 'fontsize',12,'color','w');
xlabel('X','fontsize',12);
set(gca,'fontsize',12,'yticklabel',[],...
    'tickdir','out','ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
ylim([y{i}(1) 0]);
% ylim([-4 4]);
print -djpeg b3d_xz_BFpeak.jpg


%% sqrt(TH3d_sq) to compare with Mashayek&Peltier, 
% where they have contours of density eigenfunction 
close all;
figure;
 set(gcf,'units','centimeters','paperunits','centimeters')
 set(gcf,'PaperType','A4');
 pp=[0.63 .9 15 24];
ps=[0 0 pp(3)/1.1 pp(4)/1.1];
 set(gcf,'paperposition',pp);
 set(gcf,'position',ps);

i=1;

ax1 = axes('position',[.14 .725 .8 .22]);
[cc,hh] = contourf(x{i},y{i},squeeze(sqrt(TH3d_sq{i})'),100);
set(hh,'edgecolor','none');
title(' $\sqrt{\left<b^2_{3d}\right>_y}$','fontname','times','fontsize',18,...
    'fontweight','normal','interpreter','latex');
caxis([0 0.1]);
% caxis([-1 1]);
% colormap(ax1,'default')

% hold on; contour(x{i},y{i},TH1y_cen{i}(:,:,k)',[.6,1],'k');
% colormap(ax1,(cbrewer('seq', 'PuBuGn', 100)));
% colorbar('position',[0.85,0.55 0.015 0.4]);

% text(1.5,-2.5,'$d=10, t=149$','interpreter','latex',...
% 'fontsize',15,'color','w');
text(1.5,-2.5,'$d=10, t=175$','interpreter','latex',...
'fontsize',15,'color','w');

% text(1.5,-2.5,'$d=10, t=188$','interpreter','latex',...
% 'fontsize',15,'color','w');

ylabel('Z','fontsize',12);
ylim([y{i}(1)+10-3.29 y{i}(1)+10+3.29]);
set(gca,'fontsize',12,'tickdir','out',...
    'ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;

i=2;
  ax1 = axes('position',[.14 .45 .8 .22]);
[cc,hh] = contourf(x{i},y{i},squeeze(sqrt(TH3d_sq{i})'),100);
set(hh,'edgecolor','none');
caxis([0 0.1]);
% col = colorbar('location','southoutside','position',[.09,0.4 0.19 0.015]);
% colormap(ax1,cm1)

% hold on; contour(x{i},y{i},TH1y_cen{i}(:,:,k)',[.6,1],'k');
% colormap(ax1,(cbrewer('seq', 'PuBuGn', 100)));
% colorbar('position',[0.85,0.55 0.015 0.4]);
% text(1.5,-5,'$d=5, t=154$','interpreter','latex',...
% 'fontsize',15,'color','w');
text(1.5,-5,'$d=5, t=180$','interpreter','latex',...
'fontsize',15,'color','w');

% text(1.5,-5,'$d=5, t=193$','interpreter','latex',...
% 'fontsize',15,'color','w');
ylabel('Z','fontsize',12);
ylim([y{i}(1)+5-3.29 y{i}(1)+5+3.29]);
set(gca,'fontsize',12,'tickdir','out','ytick',[-5:2:-1],...
    'ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;

i=3;

  ax1 = axes('position',[.14 .175 .8 .22]);
[cc,hh] = contourf(x{i},y{i},squeeze(sqrt(TH3d_sq{i})'),100);
set(hh,'edgecolor','none');
caxis([0 0.1]);
colormap(jet)
% hold on; contour(x{i},y{i},TH1y_cen{i}(:,:,k)',[.6,1],'k');
% colormap(ax1,(cbrewer('seq', 'PuBuGn', 100)));
% text(1.5,y{i}(1)*.8,'$\langle \omega_{x}^{2} \rangle_{y}$','interpreter','latex',...
% 'fontsize',12,'color','w');
ylabel('Z','fontsize',12);xlabel('X','fontsize',12);
col = colorbar('location','southoutside','position',[.14,0.09 0.8 0.015]);

set(gca,'fontsize',12,...
    'tickdir','out','ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;
% text(1.5,-5.8,'$d=3, t=224$','interpreter','latex',...
% 'fontsize',15,'color','w');
text(1.5,-5.8,'$d=3, t=252$','interpreter','latex',...
'fontsize',15,'color','w');

% text(1.5,-5.8,'$d=3, t=266$','interpreter','latex',...
% 'fontsize',15,'color','w');
ylim([y{i}(1) 0]);
print -djpeg rmsb3d_xz_BFpeak_onset.jpg

%%
close all;
figure;
 set(gcf,'units','centimeters','paperunits','centimeters')
 set(gcf,'PaperType','A4');
 pp=[0.63 1 20 27.5];
ps=[0 0 pp(3)/1.1 pp(4)/1.1];
 set(gcf,'paperposition',pp);
 set(gcf,'position',ps);
i=3;
k=1%:length(timestep)
ax1 = axes('position',[.08 .75 .4 .2]);
[cc,hh] = contourf(x{i},y{i},squeeze(omega_x_sq_mean(:,:,1,k))',100);
set(hh,'edgecolor','none');
title(sprintf('t=%.0f',tt{i}(k)),'fontweight','normal','fontname','times','fontsize',8);
caxis([0 max(max(max(omega_x_sq_mean(:,:,1,1:2))))*.5]);
hold on; contour(x{i},y{i},squeeze(mean(TH1y{i}(:,:,:,k),3))',5,'k');
colormap(ax1,(cbrewer('seq', 'BuPu', 100)));
% colorbar('position',[0.85,0.55 0.015 0.4]);
text(1.5,y{i}(1)*.8,'$\langle \omega_{x}^{2} \rangle_{y}$','interpreter','latex',...
'fontsize',12,'color','k');
ylabel('Z');
set(gca,'fontsize',12,'xticklabel',[],'tickdir','out',...
    'ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;
% title('$Re$=1000, $Pr$=1, $Ri$=0.16, $d$=3','fontname','times','fontsize',12,...
%     'fontweight','normal','interpreter','latex');
ylim([y{i}(1) 0]);

ax2 = axes('position',[.08 .53 .4 .2]);
[cc,hh] = contourf(x{i},y{i},squeeze(omega_z_sq_mean(:,:,1,k))',100);
set(hh,'edgecolor','none');
caxis([0 max(max(max(omega_z_sq_mean(:,:,1,1:2))))*.4]);
hold on; contour(x{i},y{i},squeeze(mean(TH1y{i}(:,:,:,k),3))',5,'k');
colormap(ax2,(cbrewer('seq', 'BuPu', 100)));
% colorbar('position',[0.85,0.55 0.015 0.4]);
text(1.5,y{i}(1)*.8,'$\langle \omega_{y}^{2} \rangle_{y}$','interpreter','latex',...
'fontsize',12,'color','k');
% xlabel('X'); 
ylabel('Z');
set(gca,'fontsize',12,'xticklabel',[],'tickdir','out',...
    'ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;
% title('$Re$=1000, $Pr$=1, $Ri$=0.16, $d$=3','fontname','times','fontsize',12,...
%     'fontweight','normal','interpreter','latex');
ylim([y{i}(1) 0])

ax3 = axes('position',[.08 .31 .4 .2]);
[cc,hh] = contourf(x{i},y{i},squeeze(omega_y_sq_mean(:,:,1,k))',100);
set(hh,'edgecolor','none');
caxis([0 max(max(max(omega_y_sq_mean(:,:,1,1:2))))*.5]);
hold on; contour(x{i},y{i},squeeze(mean(TH1y{i}(:,:,:,k),3))',5,'k');

colormap(ax3,(cbrewer('seq', 'BuPu', 100)));
% colorbar('position',[0.85,0.55 0.015 0.4]);
text(1.5,y{i}(1)*.8,'$\langle \omega_{z}^{2} \rangle_{y}$','interpreter','latex',...
'fontsize',12,'color','k');
% xlabel('X'); 
ylabel('Z');
set(gca,'fontsize',12,'xticklabel',[],'tickdir','out',...
    'ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;
ylim([y{i}(1) 0])

ax4 = axes('position',[.08 .09 .4 .2]);
[cc,hh] = contourf(x{i},y{i},squeeze(K3D{i}(:,:,k))',100);set(hh,'edgecolor','none');
% hold on; [cc,hh]=contour(x{i},y{i},tke_3D_m',10,'color','k');%set(cc,'color','k');
colormap(ax4,(cbrewer('seq', 'YlOrRd', 100)));
caxis([0 max(max(max(K3D{i}(:,:,1:2))))*.7]);
hold on; contour(x{i},y{i},squeeze(mean(TH1y{i}(:,:,:,k),3))',5,'k');
% colorbar('position',[0.85,0.12 0.015 0.4]);
text(1.5,y{i}(1)*.8,'$K_{3D}$','interpreter','latex',...
'fontsize',12,'color','k');
xlabel('X'); ylabel('Z');
set(gca,'fontsize',12,'tickdir','out','ticklength',[0.02 0.02],...
    'xminortick','on','yminortick','on');box on;
ylim([y{i}(1) 0]);

k=2;
ax1 = axes('position',[0.5 .75 .4 .2]);
[cc,hh] = contourf(x{i},y{i},squeeze(omega_x_sq_mean(:,:,k))',100);
set(hh,'edgecolor','none');
title(sprintf('t=%.0f',tt{i}(k)),'fontweight','normal','fontname','times','fontsize',8);
caxis([0 max(max(max(omega_x_sq_mean(:,:,1:2))))*.5]);
hold on; contour(x{i},y{i},squeeze(mean(TH1y{i}(:,:,:,k),3))',5,'k');
colormap(ax1,(cbrewer('seq', 'BuPu', 100)));
colorbar('position',[0.91,0.75 0.015 0.2]);
% text(1.5,y{i}(1)*.8,'$\langle \omega_{x}^{2} \rangle_{y}$','interpreter','latex',...
% 'fontsize',12,'color','w');
% ylabel('Z');
set(gca,'fontsize',12,'xticklabel',[],'yticklabel',[],...
    'tickdir','out','ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;
% title('$Re$=1000, $Pr$=1, $Ri$=0.16, $d$=3','fontname','times','fontsize',12,...
%     'fontweight','normal','interpreter','latex');
ylim([y{i}(1) 0])

ax2 = axes('position',[0.5 .53 .4 .2]);
[cc,hh] = contourf(x{i},y{i},squeeze(omega_z_sq_mean(:,:,k))',100);
set(hh,'edgecolor','none');
caxis([0 max(max(max(omega_z_sq_mean(:,:,1:2))))*.4]);
colormap(ax2,(cbrewer('seq', 'BuPu', 100)));
hold on; contour(x{i},y{i},squeeze(mean(TH1y{i}(:,:,:,k),3))',5,'k');
colorbar('position',[0.91,0.53 0.015 0.2]);
% text(1.5,y{i}(1)*.8,'$\langle \omega_{y}^{2} \rangle_{y}$','interpreter','latex',...
% 'fontsize',12,'color','w');
% xlabel('X'); 
% ylabel('Z');
set(gca,'fontsize',12,'xticklabel',[],'yticklabel',[],...
    'tickdir','out','ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;
% title('$Re$=1000, $Pr$=1, $Ri$=0.16, $d$=3','fontname','times','fontsize',12,...
%     'fontweight','normal','interpreter','latex');
ylim([y{i}(1) 0])

ax3 = axes('position',[0.5 .31 .4 .2]);
[cc,hh] = contourf(x{i},y{i},squeeze(omega_y_sq_mean(:,:,k))',100);
set(hh,'edgecolor','none');
caxis([0 max(max(max(omega_y_sq_mean(:,:,1:2))))*.5]);
hold on; contour(x{i},y{i},squeeze(mean(TH1y{i}(:,:,:,k),3))',5,'k');
colormap(ax3,(cbrewer('seq', 'BuPu', 100)));
colorbar('position',[0.91,0.31 0.015 0.2]);
% text(1.5,y{i}(1)*.8,'$\langle \omega_{z}^{2} \rangle_{y}$','interpreter','latex',...
% 'fontsize',12,'color','w');
% xlabel('X'); 
% ylabel('Z');
set(gca,'fontsize',12,'xticklabel',[],'yticklabel',[],...
    'tickdir','out','ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;
ylim([y{i}(1) 0])

ax4 = axes('position',[0.5 .09 .4 .2]);
[cc,hh] = contourf(x{i},y{i},squeeze(K3D{i}(:,:,k))',100);set(hh,'edgecolor','none');
% hold on; [cc,hh]=contour(x{i},y{i},tke_3D_m',10,'color','k');%set(cc,'color','k');
colormap(ax4,(cbrewer('seq', 'YlOrRd', 100)));
caxis([0 max(max(max(K3D{i}(:,:,1:2))))*.7]);
hold on; contour(x{i},y{i},squeeze(mean(TH1y{i}(:,:,:,k),3))',5,'k');
colorbar('position',[0.91,0.09 0.015 0.2]);

% text(1.5,y{i}(1)*.8,'$K_{3D}$','interpreter','latex',...
% 'fontsize',12,'color','w');
xlabel('X'); %ylabel('Z');
set(gca,'fontsize',12,'yticklabel',[],...
    'tickdir','out','ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
ylim([y{i}(1) 0]);

% print -djpeg omegax_K3D_0.16_d3_252_259.jpg

