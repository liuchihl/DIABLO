clear;clc;
load('stats_12_025_sta.mat','k3d','Time');

% input coefficients
% d=10, Ri=0.12
 NY=361; LY=20;
 NX=512; LX=27.93;
 NZ=128; LZ=6.98;
 Re = 1000; NU=1/Re; % Enter the Reynolds number or viscosity from input.dat
Pr=1;   kappa=NU/Pr; % Prandtl number
RI=0.12; % Enter the richardson number for each scalar
% d=10;
timestep = 400:400:12000;
addpath(genpath('/glade/u/home/liuchihl/TOOLS'));
cm_1 = load('MPL_gnuplot.rgb');

% tt = zeros(3,length(timestep));
% tke_3D_m = zeros(1,length(timestep));


base_dir=['/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/'];

filename=[base_dir 'd_2.5/Ri_0.12_0.1_2/output/'];
K3d = zeros(NX,NY,NZ,length(timestep));
TH1 = zeros(NX,NY,NZ,length(timestep));
U3d=zeros(NX,NY,NZ);V3d=zeros(NX,NY,NZ);
W3d=zeros(NX,NY,NZ);
U=zeros(NX,NY,NZ);
V=zeros(NX,NY,NZ);W=zeros(NX,NY,NZ);
for k=1:length(timestep)

  if (timestep(k)<10)
    timename=['out0000' int2str(timestep(k)) '.h5'];
  elseif (timestep(k)<100)
    timename=['out000' int2str(timestep(k)) '.h5'];
  elseif (timestep(k)<1000)
    timename=['out00' int2str(timestep(k)) '.h5'];
  elseif (timestep(k)<10000)
    timename=['out0' int2str(timestep(k)) '.h5'];
  else
    timename=['out' int2str(timestep(k)) '.h5'];
  end
varname1=['/Timestep/U'];
varname2=['/Timestep/V'];
varname3=['/Timestep/W'];
varname4=['/Timestep/TH1'];

    
if exist([filename,timename])
U=h5read([filename,timename],varname1);
V=h5read([filename,timename],varname2);
W=h5read([filename,timename],varname3);
TH1(:,:,:,k)=h5read([filename,timename],varname4);

info = h5info([filename,timename]);
%time
tt(k) = info.Groups.Attributes.Value;
NX=size(U,1); NZ=size(U,3); NY=size(U,2);
x=linspace(0,LX,NX); y=linspace(-LY/2,LY/2,NY); z=linspace(0,LZ,NZ);
% clear wx wy wz vz
% for j=1:size(U,3)
% % wx(:,:,j) = mmderiv(x,squeeze(W(:,:,j,time))) ;
% vx(:,:,j) = mmderiv(x,squeeze(V(:,:,j,time))) ;
% end
% for j=1:size(U,3)
% wy(:,:,j) = mmderiv(y,squeeze(W(:,:,j,time))')' ;
% uy(:,:,j) = mmderiv(y,squeeze(U(:,:,j,time))')' ;
% end
% for j=1:size(U,2)
% % wz(:,j,:) = mmderiv(z,squeeze(W(:,j,:,time))')' ;
% vz(:,j,:) = mmderiv(z,squeeze(V(:,j,:,time))')' ;
% end
% 
% omega_x(:,:,:,time) = wy-vz;
% % omega_y = uz-wx;
% omega_z(:,:,:,time) = vx-uy;

% buoyancy gradient
% for j=1:size(TH1,3)
% TH1y(:,:,j) = mmderiv(y,squeeze(TH1(:,:,j,time))')' ;
% end


% U_mean = mean(U,3); U_rms = sqrt((U-U_mean).^2);
% V_mean = mean(V,3); V_rms = sqrt((V-V_mean).^2);
% W_mean = mean(W,3); W_rms = sqrt((W-W_mean).^2);
end
k
U3d = U-mean(U,3);V3d = V-mean(V,3);W3d = W-mean(W,3);
K3d(:,:,:,k) = .5*(U3d.^2+V3d.^2+W3d.^2);  

end


%% plot the isosurface
% close all;
[xx,zz,yy] = meshgrid(x,z,y);
% xx = permute(xx,[2,1,3]); yy = permute(yy,[2,1,3]); 
% zz = permute(zz,[2,1,3]); 
close all;
figure('position',[50,50,800,600]);
filename=[base_dir '/movie.h5'];

plot_xy = VideoWriter('Boundary_isosurface12_25_2.avi');

plot_xy.FrameRate = 5; plot_xy.Quality = 100;
open(plot_xy);

for i=5:30
    clf;

[faces,verts,colors] = isosurface(xx,zz,yy,...
    permute(TH1(:,:,:,i),[3,1,2]),0,permute(K3d(:,:,:,i),[3,1,2]));

p1=patch('Vertices',verts,'Faces',faces,'FaceVertexCData',colors,...
    'FaceColor','interp','EdgeColor','interp');
% p1.FaceAlpha = 'flat';
view(-8,9);

% alpha(0.1)
colormap(flipud(cbrewer('div','RdYlBu',100)));
colorbar;
caxis([0,max(max(max(K3d(:,:,:,i))))*.3]);
box on; grid on;
set(gca,'fontsize',16);
xlabel('X','fontsize',18);
ylabel('Y','fontsize',18);
zlabel('Z','fontsize',18);
title(sprintf(['$Ri_0=0.12, d=2.5$, case 2, $t=%.2f$'],tt(i)),'interpreter','latex')
% ylim([-10,0])
axis([0 LX 0 LZ -10 -6]);
% xlim([0 LX]);
axes('position',[.07,.84,.15,.15]);
semilogy(Time{2},mean(k3d{2}),'linewidth',2.2);
hold on; 
plot(tt(i),mean(K3d(:,:,:,i),[1,2,3]),'.','markersize',10,'linewidth',2);
xlabel('t','fontsize',12);
ylabel('K_{3d}');
set(gca,'xtick',[0:100:500])
%   end
 drawnow
M=getframe(gcf);
 writeVideo(plot_xy,M); 

 end
% 
 close(plot_xy)

%%
%%%
% p = patch(isosurface(xx,zz,yy,permute(TH1(:,:,:,i),[3,1,2]),1));
% isonormals(xx,zz,yy,permute(TH1(:,:,:,i),[3,1,2]),p);
% p.FaceColor = 'blue';
% p.EdgeColor = 'none';
% view(25,15);
% camlight 
% lighting gouraud
% axis([0 LX 0 LZ -5 0]);
% alpha(.5);
% box on; 
% %%%
% title(sprintf('Time = %.0f',tt(i)),'fontsize',12,'fontname','times')
% view(25,43);
% axis vis3d;
% hold on;
% slice(xx,zz,yy,permute(omega_x,[3,1,2]),[],120,[]);shading flat
% slice(xx,zz,yy,permute(omega_x,[3,1,2]),0,[],[]);shading flat
%  
% % slice(xx,zz,yy,permute(omega_x,[3,1,2]),[],[],-4);shading flat
% 
% colormap(flipud(cbrewer('div', 'RdYlBu', 100)));

% colormap(flipud(cm_1));
% xlabel('X','fontsize',12);
% ylabel('Y','fontsize',12);
% zlabel('Z','fontsize',12);
% axis([0 LX 0 LZ -5 0]);
% % alpha(.9);
% box on; 
% caxis([-5 5]*1e-1);
% % co=colorbar('position',[.37 0.55 .015 .25]);
% co=colorbar;
% ylabel(co,'\omega_{x}','fontsize',12);
% 
% grid on
% print -djpeg 3D_120_sym_191_blue.jpg