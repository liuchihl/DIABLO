% This script illustrates how to load in 3D model data and calculate 
% some basic diagnostics
% Run after readmean.m
clear;clc;
% input coefficients
fname = '/LY_3/';  
LX=15.7; LY=10; LZ=3;
timestep = [18000,20000,22000,24000];
addpath(genpath('D:/NTU resource/TOOLS'));
cm_1 = load('MPL_gnuplot.rgb');

% tt = zeros(3,length(timestep));
% tke_3D_m = zeros(1,length(timestep));
for i = 1 
if     i==1    
    base_dir=['D:/Box Sync/OSU/RESEARCH/DIABLO/coarse_resolution/Optimize_LY',fname];
end
filename=[base_dir];
for k=timestep

  if (k<10)
    timename=['out0000' int2str(k) '.h5'];
  elseif (k<100)
    timename=['out000' int2str(k) '.h5'];
  elseif (k<1000)
    timename=['out00' int2str(k) '.h5'];
  elseif (k<10000)
    timename=['out0' int2str(k) '.h5'];
  else
    timename=['out' int2str(k) '.h5'];
  end
varname1=['/Timestep/U'];
varname2=['/Timestep/V'];
varname3=['/Timestep/W'];
varname4=['/Timestep/TH1'];
% index for time
if k==timestep(1)
   time = 1;
elseif k==timestep(2)
    time = 2;
elseif k==timestep(3)
    time = 3;
else k==timestep(4)
    time = 4;

end
    
if exist([filename,timename])
U(:,:,:,time)=h5read([filename,timename],varname1);
V(:,:,:,time)=h5read([filename,timename],varname2);
W(:,:,:,time)=h5read([filename,timename],varname3);
TH1(:,:,:,time)=h5read([filename,timename],varname4);

info = h5info([filename,timename]);
%time
tt(time) = info.Groups.Attributes.Value;
NX=size(U,1); NZ=size(U,3); NY=size(U,2);
x=linspace(0,LX,NX); y=linspace(-LY/2,LY/2,NY); z=linspace(0,LZ,NZ);

drhodz1=0.0;
drhodz2=0.0;
clear wx wy wz vz
for j=1:size(U,3)
% wx(:,:,j) = mmderiv(x,squeeze(W(:,:,j,time))) ;
vx(:,:,j) = mmderiv(x,squeeze(V(:,:,j,time))) ;
end
for j=1:size(U,3)
wy(:,:,j) = mmderiv(y,squeeze(W(:,:,j,time))')' ;
uy(:,:,j) = mmderiv(y,squeeze(U(:,:,j,time))')' ;
end
for j=1:size(U,2)
% wz(:,j,:) = mmderiv(z,squeeze(W(:,j,:,time))')' ;
vz(:,j,:) = mmderiv(z,squeeze(V(:,j,:,time))')' ;
end

omega_x(:,:,:,time) = wy-vz;
% omega_y = uz-wx;
omega_z(:,:,:,time) = vx-uy;


% Add the background buoyancy gradient
% for j=1:NZ
%   TH1(:,:,j)=TH1(:,:,j)+drhodz1*z(j);
% %  TH2(:,:,k)=TH2(:,:,k)+drhodz2*z(k);
% end

% calculate tke 3D
% U_mean = mean(U,3); U_rms = sqrt((U-U_mean).^2);
% V_mean = mean(V,3); V_rms = sqrt((V-V_mean).^2);
% W_mean = mean(W,3); W_rms = sqrt((W-W_mean).^2);
% tke_3D = .5*(U_rms.^2+V_rms.^2+W_rms.^2);
% tke_3D_m(i,k/inv) = mean(mean(mean(tke_3D)));
end

end
end
%% plot the isosurface
close all;
[xx,yy,zz] = meshgrid(x,y,z);
xx = permute(xx,[2,1,3]); yy = permute(yy,[2,1,3]); 
zz = permute(zz,[2,1,3]); 
figure('position',[50,50,800,600]);

axes('position',[0.02 0.55 .35 .35]);
i= 1;
[faces,verts,colors] = isosurface(xx,zz,yy,TH1(:,:,:,i),1,omega_x(:,:,:,i));
patch('Vertices',verts,'Faces',faces,'FaceVertexCData',colors,...
    'FaceColor','interp','EdgeColor','interp');
hold on;
title(sprintf('Time = %.0f',tt(i)),'fontsize',12,'fontname','times')
view(25,43)
axis vis3d;
colormap(flipud(cm_1))
xlabel('X','fontsize',12);
ylabel('Y','fontsize',12);
zlabel('Z','fontsize',12);
axis([0 LX 0 LZ -4 0]);
alpha(.5)
box on; 
caxis([-.5 .5]*1e-3);
co=colorbar('position',[.37 0.55 .015 .25]);
ylabel(co,'\omega_{x}','fontsize',12);

axes('position',[0.5 0.55 .35 .35]);
i= 2;
[faces,verts,colors] = isosurface(xx,zz,yy,TH1(:,:,:,i),1,omega_x(:,:,:,i));
patch('Vertices',verts,'Faces',faces,'FaceVertexCData',colors,...
    'FaceColor','interp','EdgeColor','interp');
title(sprintf('Time = %.0f',tt(i)),'fontsize',12,'fontname','times')
view(25,43)
axis vis3d;
colormap(flipud(cm_1))
xlabel('X','fontsize',12);
ylabel('Y','fontsize',12);
zlabel('Z','fontsize',12);
axis([0 LX 0 LZ -4 0]);
alpha(.5)
box on; 
caxis([-2 2]*1e-3);
co=colorbar('position',[.85 0.55 .015 .25]);
ylabel(co,'\omega_{x}','fontsize',12);

axes('position',[0.02 0.1 .35 .35]);
i= 3;
[faces,verts,colors] = isosurface(xx,zz,yy,TH1(:,:,:,i),1,omega_x(:,:,:,i));
patch('Vertices',verts,'Faces',faces,'FaceVertexCData',colors,...
    'FaceColor','interp','EdgeColor','interp');
title(sprintf('Time = %.0f',tt(i)),'fontsize',12,'fontname','times')
view(25,43)
axis vis3d;
colormap(flipud(cm_1))
xlabel('X','fontsize',12);
ylabel('Y','fontsize',12);
zlabel('Z','fontsize',12);
axis([0 LX 0 LZ -4 0]);
alpha(.5)
box on; 
caxis([-2 2]*1e-2);
co=colorbar('position',[.37 0.1 .015 .25]);
ylabel(co,'\omega_{x}','fontsize',12);

axes('position',[0.5 0.1 .35 .35]);
i= 4;
[faces,verts,colors] = isosurface(xx,zz,yy,TH1(:,:,:,i),1,omega_x(:,:,:,i));
patch('Vertices',verts,'Faces',faces,'FaceVertexCData',colors,...
    'FaceColor','interp','EdgeColor','interp');
title(sprintf('Time = %.0f',tt(i)),'fontsize',12,'fontname','times')
view(25,43)
axis vis3d;
colormap((cm_1));
% cmocean('deep')
xlabel('X','fontsize',12);
ylabel('Y','fontsize',12);
zlabel('Z','fontsize',12);
axis([0 LX 0 LZ -4 0]);
alpha(.8)
box on; 
co=colorbar('position',[.85 0.1 .015 .25]);
ylabel(co,'\omega_{x}','fontsize',12);
caxis([-2 2]*1e-2);

% print -djpeg 3D_3.jpg
% axis equal
% colormap hot



%% 

% clear; 
% 
% [x,y,z,v] = flow; 
% [faces,verts,colors] = isosurface(x,y,z,v,-3,v); 
% patch('Vertices',verts,'Faces',faces,'FaceVertexCData',colors,...
%     'FaceColor','interp','EdgeColor','interp')
% view(30,-15)
% axis vis3d
% colormap copper
% 
% % alpha(.2)
% 
% caxis([-3.1 -2.8])

%%





%%
% figure; 
% for i=1
% plot(tt(i,:),tke_3D_m(i,:),'o-');
% hold on;
% end
% le = legend('Ly = 4','Ly = 5','Ly = 6','Ly = 7','Ly = 8');
% le = legend('Ly = 4','Ly = 6','Ly = 7');
% set(le,'fontsize',12);
% xlabel('time','fontsize',12); ylabel('K_3_D','fontsize',12);
% % Calculate the x-average velocity
% ume=squeeze(mean(U,1));
% vme=squeeze(mean(V,1));
% wme=squeeze(mean(W,1));
% thme1=squeeze(mean(TH1,1));
% %thme2=squeeze(mean(TH2,1));
% 
% % Calculate correlation terms
% for k=1:NZ
%   uw_mean(:,k)=mean((U(:,:,k)-U_BT(k)).*(W(:,:,k)-W_BT(k)),1);
%   uu_mean(:,k)=mean((U(:,:,k)-U_BT(k)).*(U(:,:,k)-U_BT(k)),1);
%   ww_mean(:,k)=mean((W(:,:,k)-W_BT(k)).*(W(:,:,k)-W_BT(k)),1);
%   uw_BT(:,k)=trapz(gyf,uw_mean(:,k))/(gyf(end)-gyf(1));
%   uu_BT(:,k)=trapz(gyf,uu_mean(:,k))/(gyf(end)-gyf(1));
%   ww_BT(:,k)=trapz(gyf,ww_mean(:,k))/(gyf(end)-gyf(1));
% end
% 
% % calculate the mean shear
% dudy=zeros(size(ume));
% dwdy=zeros(size(wme));
% for j=2:NY-1
%     dudy(j,:)=(ume(j+1,:)-ume(j-1,:))/(gyf(j+1)-gyf(j-1));
%     dwdy(j,:)=(wme(j+1,:)-wme(j-1,:))/(gyf(j+1)-gyf(j-1));
% end
% 
% % Calculate the local shear
% dUdy=zeros(size(U));
% for j=2:NY-1
% dUdy(:,j,:)=(U(:,j+1,:)-U(:,j-1,:))/(gyf(j+1)-gyf(j-1));
% end


