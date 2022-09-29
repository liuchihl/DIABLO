% This script calculates the K3d budget and with plots in the x-z plane. 
% We can see how each terms distribute and further understand the mechanism 

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

timestep = [28000;29000;38000]; % the timestep where turbulence has developed initially

%    2 peaks of APE:
% 	d=10: t=153 (22000), 200 (29000)
% 	d=5: t=154 (23000), 201 (30000)
% 	d=3: t=221 (32000), 277 (40000)
% timestep = [36000];
addpath(genpath('/glade/u/home/liuchihl/TOOLS'));

for i = 1:3
U=zeros(NX,NY(i),NZ);V=zeros(NX,NY(i),NZ);
W=zeros(NX,NY(i),NZ);TH1=zeros(NX,NY(i),NZ);

U3d=zeros(NX,NY(i),NZ);V3d=zeros(NX,NY(i),NZ);W3d=zeros(NX,NY(i),NZ);
TH3d=zeros(NX,NY(i),NZ);
U2d=zeros(NX,NY(i));V2d=zeros(NX,NY(i));W2d=zeros(NX,NY(i));

 U3d_x=zeros(NX,NY(i),NZ);U3d_y=zeros(NX,NY(i),NZ);U3d_z=zeros(NX,NY(i),NZ);
 V3d_y=zeros(NX,NY(i),NZ);V3d_x=zeros(NX,NY(i),NZ);V3d_z=zeros(NX,NY(i),NZ);
 W3d_x=zeros(NX,NY(i),NZ);W3d_y=zeros(NX,NY(i),NZ);W3d_z=zeros(NX,NY(i),NZ);
 
 SP_2d{i}=zeros(NX,NY(i));SP_background{i}=zeros(NX,NY(i));
 BF3d{i}=zeros(NX,NY(i));D3d{i}=zeros(NX,NY(i)); SP_sheardeform{i}=zeros(NX,NY(i));

K3d{i}=zeros(NX,NY(i));


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
tt{i} = info.Groups.Attributes.Value;
% NX=size(U,1); NZ=size(U,3); NY=size(U,2);
x{i}=linspace(0,LX(i),NX); y{i}=linspace(-LY(i)/2,LY(i)/2,NY(i)); 
z{i}=linspace(0,LZ(i),NZ);

% calculate 3D perturbation kinetic energy (Smyth et al., 2005)
U3d = U-mean(U,3);V3d = V-mean(V,3);W3d = W-mean(W,3);
TH3d = TH1-mean(TH1,3);
K3d{i} = mean(.5*(U3d.^2+V3d.^2+W3d.^2),3);  %

for j=1:NY(i)
% TH3d_x(:,j,:)= gradient(squeeze(TH3d(:,j,:))',mean(diff(x{i})))';
 U3d_x(:,j,:) = gradient(squeeze(U3d(:,j,:))',mean(diff(x{i})))';
 V3d_x(:,j,:) = gradient(squeeze(V3d(:,j,:))',mean(diff(x{i})))';
 W3d_x(:,j,:) = gradient(squeeze(W3d(:,j,:))',mean(diff(x{i})))';
end
for j=1:NX
 V3d_y(j,:,:) = gradient(squeeze(V3d(j,:,:))',mean(diff(y{i})))';
 U3d_y(j,:,:) = gradient(squeeze(U3d(j,:,:))',mean(diff(y{i})))';
 W3d_y(j,:,:) = gradient(squeeze(W3d(j,:,:))',mean(diff(y{i})))';
 U3d_z(j,:,:) = gradient(squeeze(U3d(j,:,:)),mean(diff(z{i})));
 V3d_z(j,:,:) = gradient(squeeze(V3d(j,:,:)),mean(diff(z{i})));
 W3d_z(j,:,:) = gradient(squeeze(W3d(j,:,:)),mean(diff(z{i})));
end

end

% calculate 2D perturbation kinetic energy (Smyth et al., 2005)
 U2d = mean(U,3)-mean(mean(U,1),3);
 V2d = mean(V,3)-mean(mean(V,1),3);
 W2d = mean(W,3)-mean(mean(W,1),3);
% K2d{i}(k) = mean(mean(mean(.5*(U2d.^2+V2d.^2+W2d.^2),3)));  %
%
% % % calculate the K3d energy budget
% % shear production: extract energy from the mean background field
% %
 SP_background{i} = -mean(...
     U3d.*V3d.*gradient(mean(mean(U,1),3),mean(diff(y{i}))) ,3 );

 % shear production: extract energy from the 2D KH flow field
 SP_2d{i} = -mean(...
     U3d.^2.*gradient(U2d',mean(diff(x{i})))'+...
     V3d.^2.*gradient(V2d,mean(diff(y{i}))) ,3);
 % shearing deformation
 SP_sheardeform{i} = -mean(...
     U3d.*V3d.*gradient(U2d,mean(diff(y{i})))+gradient(V2d',mean(diff(x{i})))' ,3);

 % buoyancy flux
 BF3d{i} = RI*mean(...
     TH3d.*V3d ,3);
 % dissipation
 D3d{i} = -1/Re*mean(...
      2*(U3d_x.^2+V3d_y.^2+W3d_z.^2)+U3d_y.^2+U3d_z.^2+...
      V3d_x.^2+V3d_z.^2+W3d_x.^2+W3d_y.^2+...
      2*(U3d_z.*W3d_x+U3d_y.*V3d_x+V3d_z.*W3d_y) ,3);
k

i
end

end
%tke_3D_m{i}
%SP_background{i}+SP_2d{i}+SP_sheardeform{i}+BF3d{i}+D3d{i}
 save K3D_budget_xy SP_background SP_2d SP_sheardeform BF3d D3d K3d x y z tt






%% %% plot the x-z of the k3d budget  TH3d^2/2
cm_1 = load('MPL_gnuplot.rgb');
cm_eps = load('NCV_rainbow2.rgb');

% load K3D_budget_xy.mat
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
caxis([0 0.05]);
% caxis([-1 1]);
 colormap(ax1,(cm_eps/255));
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
caxis([-8e-3 .015]);
colormap(ax2,jet);
title('$SP_{background}$','interpreter','latex',...
'fontsize',18,'color','k');
set(gca,'fontsize',12,'yticklabel',[],'tickdir','out',...
    'ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;
ylim([y{i}(1)+10-3.29 y{i}(1)+10+3.29]);


ax3 = axes('position',[.62 .725 .26 .22]);
  [cc,hh] = contourf(x{i},y{i},squeeze(D3d{i})',100);set(hh,'edgecolor','none');
% colormap(ax3,flipud((cbrewer('seq', 'YlGn', 100))));
 colormap(ax3,flipud(parula));
caxis([-2.5,0]*1e-3);
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
caxis([0 0.05]);
% caxis([-1 1]);
% col = colorbar('location','southoutside','position',[.08,0.4 0.26 0.015]);
 colormap(ax1,(cm_eps/255));
text(1.5,-5,'$d=5, t=193$','interpreter','latex',...
'fontsize',15,'color','w');
ylabel('Z','fontsize',12);
ylim([y{i}(1)+5-3.29 y{i}(1)+5+3.29]);
set(gca,'fontsize',12,'tickdir','out',...
    'ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;


ax2 = axes('position',[.35 .45 .26 .22]);
    [cc,hh] = contourf(x{i},y{i},squeeze(SP_background{i})',100);
set(hh,'edgecolor','none');
caxis([-8e-3 .015]);
colormap(ax2,jet);
set(gca,'fontsize',12,'yticklabel',[],'tickdir','out',...
    'ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;
ylim([y{i}(1)+5-3.29 y{i}(1)+5+3.29]);


ax3 = axes('position',[.62 .45 .26 .22]);
  [cc,hh] = contourf(x{i},y{i},squeeze(D3d{i})',100);set(hh,'edgecolor','none');
% colormap(ax3,flipud((cbrewer('seq', 'YlGn', 100))));
 colormap(ax3,flipud(parula));

caxis([-2.5,0]*1e-3);
set(gca,'fontsize',12,'tickdir','out','ticklength',[0.02 0.02],...
    'xminortick','on','yminortick','on');box on;
ylim([y{i}(1)+5-3.29 y{i}(1)+5+3.29]);
set(gca,'fontsize',12,'yticklabel',[],'tickdir','out',...
    'ticklength',[0.02 0.02],'xminortick','on','yminortick','on');

i=3;
 ax1 = axes('position',[.08 .175 .26 .22]);
[cc,hh] = contourf(x{i},y{i},squeeze(K3d{i}'),100);
set(hh,'edgecolor','none');
caxis([0 0.05]);
 colormap(ax1,(cm_eps/255));
ylabel('Z','fontsize',12);xlabel('X','fontsize',12);
col = colorbar('location','southoutside','position',[.1,0.06 0.22 0.015]);
set(gca,'fontsize',12,...
    'tickdir','out','ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;
text(1.5,-5.8,'$d=3, t=266$','interpreter','latex',...
'fontsize',15,'color','w');
ylim([y{i}(1) 0])
 
ax2 = axes('position',[.35 .175 .26 .22]);
[cc,hh] = contourf(x{i},y{i},squeeze(SP_background{i})',100);
set(hh,'edgecolor','none');
col = colorbar('location','southoutside','position',[.37,0.06 0.22 0.015]);
caxis([-8e-3 .015]);
colormap(ax2,jet);
xlabel('X','fontsize',12);
set(gca,'fontsize',12,'yticklabel',[],...
    'tickdir','out','ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;
ylim([y{i}(1) 0]);


ax3 = axes('position',[.62 .175 .26 .22]);
[cc,hh] = contourf(x{i},y{i},squeeze(D3d{i})',100);set(hh,'edgecolor','none');
caxis([-2.5,0]*1e-3);
% colormap(ax3,flipud((cbrewer('seq', 'YlGn', 100))));
 colormap(ax3,flipud(parula));
col = colorbar('location','southoutside','position',[.64,0.06 0.22 0.015]);
xlabel('X','fontsize',12);
set(gca,'fontsize',12,'yticklabel',[],...
    'tickdir','out','ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
ylim([y{i}(1) 0]);
print -djpeg K3d_SPshear_eps.jpg

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
caxis([-0.01 0.01]);
% caxis([-1 1]);
 colormap(ax1,(cm_eps/255));
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
caxis([-.3 .3]);
colormap(ax2,jet);
title('$SP_{shear}$','interpreter','latex',...
'fontsize',18,'color','k');
set(gca,'fontsize',12,'yticklabel',[],'tickdir','out',...
    'ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;
ylim([y{i}(1)+10-3.29 y{i}(1)+10+3.29]);


ax3 = axes('position',[.62 .725 .26 .22]);
  [cc,hh] = contourf(x{i},y{i},squeeze(BF3d{i})',100);set(hh,'edgecolor','none');
% colormap(ax3,flipud((cbrewer('seq', 'YlGn', 100))));
 colormap(ax3,flipud(parula));
caxis([-2,2]*1e-3);
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
caxis([-0.01 0.01]);
% caxis([-1 1]);
% col = colorbar('location','southoutside','position',[.08,0.4 0.26 0.015]);
 colormap(ax1,(cm_eps/255));
text(1.5,-5,'$d=5, t=193$','interpreter','latex',...
'fontsize',15,'color','w');
ylabel('Z','fontsize',12);
ylim([y{i}(1)+5-3.29 y{i}(1)+5+3.29]);
set(gca,'fontsize',12,'tickdir','out',...
    'ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;


ax2 = axes('position',[.35 .45 .26 .22]);
    [cc,hh] = contourf(x{i},y{i},squeeze(SP_sheardeform{i})',100);
set(hh,'edgecolor','none');
caxis([-.3 .3]);
colormap(ax2,jet);
set(gca,'fontsize',12,'yticklabel',[],'tickdir','out',...
    'ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;
ylim([y{i}(1)+5-3.29 y{i}(1)+5+3.29]);


ax3 = axes('position',[.62 .45 .26 .22]);
  [cc,hh] = contourf(x{i},y{i},squeeze(BF3d{i})',100);set(hh,'edgecolor','none');
% colormap(ax3,flipud((cbrewer('seq', 'YlGn', 100))));
 colormap(ax3,flipud(parula));

caxis([-2,2]*1e-3);
set(gca,'fontsize',12,'tickdir','out','ticklength',[0.02 0.02],...
    'xminortick','on','yminortick','on');box on;
ylim([y{i}(1)+5-3.29 y{i}(1)+5+3.29]);
set(gca,'fontsize',12,'yticklabel',[],'tickdir','out',...
    'ticklength',[0.02 0.02],'xminortick','on','yminortick','on');

i=3;
 ax1 = axes('position',[.08 .175 .26 .22]);
[cc,hh] = contourf(x{i},y{i},squeeze(SP_2d{i}'),100);
set(hh,'edgecolor','none');
caxis([-0.01 0.01]);
 colormap(ax1,(cm_eps/255));
ylabel('Z','fontsize',12);xlabel('X','fontsize',12);
col = colorbar('location','southoutside','position',[.1,0.06 0.22 0.015]);
set(gca,'fontsize',12,...
    'tickdir','out','ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;
text(1.5,-5.8,'$d=3, t=266$','interpreter','latex',...
'fontsize',15,'color','w');
ylim([y{i}(1) 0])
 
ax2 = axes('position',[.35 .175 .26 .22]);
[cc,hh] = contourf(x{i},y{i},squeeze(SP_sheardeform{i})',100);
set(hh,'edgecolor','none');
col = colorbar('location','southoutside','position',[.37,0.06 0.22 0.015]);
caxis([-.3 .3]);
set(col,'xtick',[-.3:0.3:.3],'xticklabel',[-.3:.3:.3])
colormap(ax2,jet);
xlabel('X','fontsize',12);
set(gca,'fontsize',12,'yticklabel',[],...
    'tickdir','out','ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
box on;
ylim([y{i}(1) 0]);


ax3 = axes('position',[.62 .175 .26 .22]);
[cc,hh] = contourf(x{i},y{i},squeeze(BF3d{i})',100);set(hh,'edgecolor','none');
caxis([-2,2]*1e-3);
% colormap(ax3,flipud((cbrewer('seq', 'YlGn', 100))));
 colormap(ax3,flipud(parula));
col = colorbar('location','southoutside','position',[.64,0.06 0.22 0.015]);
xlabel('X','fontsize',12);
set(gca,'fontsize',12,'yticklabel',[],...
    'tickdir','out','ticklength',[0.02 0.02],'xminortick','on','yminortick','on');
ylim([y{i}(1) 0]);
print -djpeg SPkh_SPshear_BF.jpg



