% This script illustrates how to load in 3D model data and calculate 
% some basic diagnostics
% Run after readmean.m
%clear;%clc; 

% input coefficients
fname1 = 'd_10/Ri_0.16/'; 
% fname2= 'd_5/Ri_0.16/';  
%fname3 = 'd_3/Ri_0.16/';  
% fname1= 'd_3/Ri_0.12_0.01_20/';
% In DIABLO coordinate
%NY=[721,541,397];         % Here, NY should match the value in grid_def.all
%LY = [19.86,14.75,13.16];     %vertical length
%NX=1024;
%LX=[28.2 27.9 28.8];
%LZ=[7.05,6.98,7.2];  
%NZ=256;
RI=0.16; Re=1000;Pr=1;

%%% d=3, LY=20
%Ri=0.12
LX=27.76; NX=512;
LY=20;    NY=361; 	 
LZ=6.94;  NZ=128;
timestep=[500:500:20000];

% %Ri=0.16
%LX=28.61; NX=1024;
%LY=20;   NY=712;
%LZ=7.15;  NZ=256;

%%% d=10, LY=20
% %Ri=0.12
%LX=28.28; NX=1024;
%LY=20;   NY=721;
%LZ=7.07;  NZ=256;
% %Ri=0.16
%LX=27.73; NX=1024;
%LY=20;   NY=739;
%LZ=6.93;  NZ=256;

% test (small)
%fname1= 'KH_test4';
%LX=30; LY=20; LZ=4;
%NX=128;NY=65; NZ=16;
%timestep=[200:200:1000];

%timestep = [36000];
%   timestep = [1000:1000:60000];
%timestep=[200:200:24000];
% timestep = [36000];
addpath(genpath('/glade/u/home/liuchihl/TOOLS'));
cm_1 = load('MPL_gnuplot.rgb');
for i = 1
U=zeros(NX,NY(i),NZ);V=zeros(NX,NY(i),NZ);
W=zeros(NX,NY(i),NZ);TH1=zeros(NX,NY(i),NZ);

U3d=zeros(NX,NY(i),NZ);V3d=zeros(NX,NY(i),NZ);W3d=zeros(NX,NY(i),NZ);
TH3d=zeros(NX,NY(i),NZ);
U2d=zeros(NX,NY(i));V2d=zeros(NX,NY(i));W2d=zeros(NX,NY(i));

 U3d_x=zeros(NX,NY(i),NZ);U3d_y=zeros(NX,NY(i),NZ);U3d_z=zeros(NX,NY(i),NZ);
 V3d_y=zeros(NX,NY(i),NZ);V3d_x=zeros(NX,NY(i),NZ);V3d_z=zeros(NX,NY(i),NZ);
 W3d_x=zeros(NX,NY(i),NZ);W3d_y=zeros(NX,NY(i),NZ);W3d_z=zeros(NX,NY(i),NZ);
TH3d_x=zeros(NX,NY(i),NZ); TH3d_y=zeros(NX,NY(i),NZ);
TH3d_z=zeros(NX,NY(i),NZ);

 SP_background{i}=zeros(size(timestep));
 SP_2d{i}=zeros(size(timestep));
 BF3d{i}=zeros(size(timestep)); D3d{i}=zeros(size(timestep)); 
 SP_sheardeform{i}=zeros(size(timestep));
K3d{i}=zeros(size(timestep)); K2d{i}=zeros(size(timestep));
dTHzdy{i}=zeros(NX,NY(i)); dTHzdx{i}=zeros(NX,NY(i));
LHS{i}=zeros(size(timestep));
D3d_p{i}=zeros(size(timestep));
%if     i==1    
%    base_dir=['/glade/campaign/univ/uosc0024/diablo_2.0/Large_case/',fname1,'output/'];
%    %base_dir=['/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_10/Ri_0.16_0.5/','output/'];
%elseif i==2
%    base_dir=['/glade/campaign/univ/uosc0024/diablo_2.0/Large_case/',fname2,'output/'];
%else
%    %base_dir=['/glade/campaign/univ/uosc0024/diablo_2.0/Large_case/',fname3,'output/'];
%    base_dir=['/glade/campaign/univ/uosc0024/diablo_2.0/Large_case/',fname3,'output/'];
%
%end
%base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/KH_test4/output/';
 base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_10/butterfly/Ri_0.16_0.05_1/output/';
filename=[base_dir];
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
U(:,:,:)=h5read([filename,timename],varname1);
V(:,:,:)=h5read([filename,timename],varname2);
W(:,:,:)=h5read([filename,timename],varname3);
TH1(:,:,:)=h5read([filename,timename],varname4);

info = h5info([filename,timename]);
%time
tt{i}(k) = info.Groups.Attributes.Value;
% NX=size(U,1); NZ=size(U,3); NY=size(U,2);
x=linspace(0,LX(i),NX); y=gyf;%linspace(-LY(i)/2,LY(i)/2,NY(i)); 
z=linspace(0,LZ(i),NZ);


% vorticities
% for j=1:size(U,3)
% % wx(:,:,j) = mmderiv(x,squeeze(W(:,:,j,time))) ;
% vx(:,:,j) = mmderiv(x,squeeze(V(:,:,j,time))) ;
% end
% for j=1:size(U,3)
% wy(:,:,j) = mmderiv(y,squeeze(W(:,:,j,k))')' ;
% uy(:,:,j) = mmderiv(y,squeeze(U(:,:,j,k))')' ;
% end
% for j=1:size(U,2)
% % wz(:,j,:) = mmderiv(z,squeeze(W(:,j,:,time))')' ;
% vz(:,j,:) = mmderiv(z,squeeze(V(:,j,:,k))')' ;
% end
% 
% omega_x(:,:,:,k) = wy-vz;
% omega_y = uz-wx;
% omega_z(:,:,:,time) = vx-uy;

% % buoyancy gradient
% for j=1:size(TH1,3)
% TH1y(:,:,j) = mmderiv(y,squeeze(TH1(:,:,j,time))')' ;
% end

% Add the background buoyancy gradient
% for j=1:NZ
%   TH1(:,:,j)=TH1(:,:,j)+drhodz1*z(j);
% %  TH2(:,:,k)=TH2(:,:,k)+drhodz2*z(k);
% end


% U_mean = mean(U,3); U_rms = sqrt((U-U_mean).^2);
% V_mean = mean(V,3); V_rms = sqrt((V-V_mean).^2);
% W_mean = mean(W,3); W_rms = sqrt((W-W_mean).^2);
end

% calculate 3D perturbation kinetic energy (Smyth et al., 2005)
U3d = U-mean(U,3);V3d = V-mean(V,3);W3d = W-mean(W,3);
TH3d = TH1-mean(TH1,3);
K3d{i}(k) = mean(mean(mean(.5*(U3d.^2+V3d.^2+W3d.^2))));  %
%K3d_umez{i}(k) = mean(mean(mean(.5*((U-UME_Z).^2+().^2+W3d.^2))));
% 3D perturbation derivatives

%for j=1:NY(i)
% U3d_x(:,j,:) = mmderiv(x,squeeze(U3d(:,j,:)));
% V3d_x(:,j,:) = mmderiv(x,squeeze(V3d(:,j,:)));
% W3d_x(:,j,:) = mmderiv(x,squeeze(W3d(:,j,:)));
% U3d_x(:,j,:) = gradient(squeeze(U3d(:,j,:))',mean(diff(x)))';
% V3d_x(:,j,:) = gradient(squeeze(V3d(:,j,:))',mean(diff(x)))';
% W3d_x(:,j,:) = gradient(squeeze(W3d(:,j,:))',mean(diff(x)))';
%TH3d_x(:,j,:)= gradient(squeeze(TH3d(:,j,:))',mean(diff(x)))';
%end

%for j=1:NX
% V3d_y(j,:,:) = mmderiv(y,squeeze(V3d(j,:,:)));
% U3d_y(j,:,:) = mmderiv(y,squeeze(U3d(j,:,:)));
% W3d_y(j,:,:) = mmderiv(y,squeeze(W3d(j,:,:)));
% V3d_y(j,:,:) = gradient(squeeze(V3d(j,:,:))',mean(diff(y)))';
% U3d_y(j,:,:) = gradient(squeeze(U3d(j,:,:))',mean(diff(y)))';
% W3d_y(j,:,:) = gradient(squeeze(W3d(j,:,:))',mean(diff(y)))';
%TH3d_y(j,:,:) = gradient(squeeze(TH3d(j,:,:))',mean(diff(y)))';

% U3d_z(j,:,:) = mmderiv(z,squeeze(U3d(j,:,:))')';
% V3d_z(j,:,:) = mmderiv(z,squeeze(V3d(j,:,:))')';
% W3d_z(j,:,:) = mmderiv(z,squeeze(W3d(j,:,:))')';
% U3d_z(j,:,:) = gradient(squeeze(U3d(j,:,:)),mean(diff(z)));
% V3d_z(j,:,:) = gradient(squeeze(V3d(j,:,:)),mean(diff(z)));
% W3d_z(j,:,:) = gradient(squeeze(W3d(j,:,:)),mean(diff(z)));
%TH3d_z(j,:,:) = gradient(squeeze(TH3d(j,:,:)),mean(diff(z)));

%end


% calculate 2D perturbation kinetic energy (Smyth et al., 2005)
 U2d = mean(U,3)-mean(mean(U,1),3);
 V2d = mean(V,3)-mean(mean(V,1),3);
 W2d = mean(W,3)-mean(mean(W,1),3);
 K2d{i}(k) = mean(mean(mean(.5*(U2d.^2+V2d.^2+W2d.^2),3)));  %
%
% % % calculate the K3d energy budget
% % shear production: extract energy from the mean background field
%
% SP_background{i}(k) = -mean(mean(mean(...
%     U3d.*V3d.*repmat(gradient(mean(mean(U,1),3),mean(diff(y))),NX,1,NZ) ))); 
 		
 % shear production: extract energy from the 2D KH flow field
% SP_2d{i}(k) = -mean(mean(mean(...
%     U3d.^2.*repmat(gradient(U2d',mean(diff(x)))',1,1,NZ)+...
%     V3d.^2.*repmat(gradient(V2d,mean(diff(y))),1,1,NZ) )));
 % shearing deformation
% SP_sheardeform{i}(k) = -mean(mean(mean(...
%     U3d.*V3d.*(repmat(gradient(U2d,mean(diff(y))),1,1,NZ)+repmat(gradient(V2d',mean(diff(x)))',1,1,NZ) ) )));
		
 % buoyancy flux
% BF3d{i}(k) = RI*mean(mean(mean(...
%     TH3d.*V3d)));
 % dissipation
% D3d{i}(k) = -1/Re*mean(mean(mean(...
%      2*(U3d_x.^2+V3d_y.^2+W3d_z.^2)+U3d_y.^2+U3d_z.^2+...
%      V3d_x.^2+V3d_z.^2+W3d_x.^2+W3d_y.^2+...
%      2*(U3d_z.*W3d_x+U3d_y.*V3d_x+V3d_z.*W3d_y) )));
k

% calculate the buoyancy variance budget equation
% dissipation rate (epsilon_p) 
%dTHzdy{i} = gradient(squeeze(mean(TH1,3)),mean(diff(y)));  
%dTHzdx{i} = gradient(squeeze(mean(TH1,3))',mean(diff(x)))';   

%D3d_p{i}(k) = -1/(Re*Pr)*mean(mean(mean(...
%     (TH3d_x.^2+TH3d_y.^2+TH3d_z.^2) )));
% buoyancy flux
% BF3d_v{i}(k) = -mean(mean(mean(...
%     dTHzdy{i}.*TH3d.*V3d )));
% BF3d_h{i}(k) = -mean(mean(mean(...
%     dTHzdx{i}.*TH3d.*U3d )));
%LHS{i}(k) = mean(mean(mean(...
%    1/2*TH3d.^2 )));
end
% if i==1
%  save K3D_budget_d3_20 SP_background SP_2d SP_sheardeform BF3d D3d K3d tt %D3d_p TH3d dTHydy tt
%save B3d_budget_d10 D3d_p BF3d_v BF3d_h LHS tt

% elseif i==2
%  save K3D_budget_d5 SP_background SP_2d SP_sheardeform BF3d D3d K3d tt %D3d_p TH3d dTHydy tt
%save B3d_budget_d5 D3d_p BF3d_v BF3d_h LHS tt

% else
% save K3D_budget_d3_20 SP_background SP_2d SP_sheardeform BF3d D3d K3d K2d tt %D3d_p TH3d dTHydy tt
%save B3d_budget_d3 D3d_p BF3d_v BF3d_h LHS tt
%end

end
%tke_3D_m{i}
%SP_background{i}+SP_2d{i}+SP_sheardeform{i}+BF3d{i}+D3d{i}
%save VT5 VT tt
%save DT5 DT tt
%save B3d_budget_dot D3d_p BF3d LHS tt
% save K3D_budget_d10 SP_background SP_2d SP_sheardeform BF3d D3d K3d tt %D3d_p TH3d dTHydy tt
% save('omega_K3D_16_3_20000_30000','omega_x','K3D','tt','-v7.3')
%% plot the budget
%sigma=1./(2*tke_3D_m{i}).*mmderiv(tt{i},tke_3D_m{i});
%figure; 
%plot(tt{i},sigma,'linewidth',1.3); hold on;
%plot(tt{i},1./(2*tke_3D_m{i}).*SP_background{i},'linewidth',1.3);
%plot(tt{i},1./(2*tke_3D_m{i}).*SP_2d{i},'linewidth',1.3);
%plot(tt{i},1./(2*tke_3D_m{i}).*SD3d{i},'linewidth',1.3);
%plot(tt{i},1./(2*tke_3D_m{i}).*BF3d{i},'linewidth',1.3);
%plot(tt{i},1./(2*tke_3D_m{i}).*D3d{i},'linewidth',1.3);
%set(gca,'fontsize',12,'xminortick','on','yminortick','on');
%legend('\sigma_{3d}','SP_{background}','SP_{kh}','Stretching','BF','D_{3d}');
%xlim([0 tt{i}(end)])
%xlabel('time','fontsize',12); ylabel('\sigma_{3d}','fontsize',12);
% print -djpeg K3D_budget.jpg

%% plot K2D and K3D
%color=[84, 68, 187;213, 21, 21;255, 128, 0]/255;

%figure; 
%for i=1:3
%a(i)=plot(tt{i}(1:56),tke_2D_m{i}(1:56),'linewidth',1.5,'color',color(i,:)); hold on;
%b(i)=plot(tt{i}(1:56),tke_3D_m{i}(1:56),':','linewidth',1.5,'color',color(i,:));
%end
%xlabel('time');
%legend([a(1),a(2),a(3)],'d=10','d=5','d=3');
%text(150,5.5e-3,'K_{2D}'); text(200,5.5e-3,'K_{3D}');
%set(gca,'fontsize',12,'xminortick','on','yminortick','on'); 
%print -djpeg KE.jpg
%% plot the isosurface in terms of omega_x
% close all;

%[xx,zz,yy] = meshgrid(x,z,y);
% xx = permute(xx,[2,1,3]); yy = permute(yy,[2,1,3]); 
% zz = permute(zz,[2,1,3]); 
%figure('position',[50,50,800,600]);

% axes('position',[0.1 0.55 .8 .4]);
%i= 1;
%[faces,verts,colors] = isosurface(xx,zz,yy,...
%    permute(TH1(:,:,:,i),[3,1,2]),.8,permute(omega_x(:,:,:,i),[3,1,2]));
%patch('Vertices',verts,'Faces',faces,'FaceVertexCData',colors,...
%    'FaceColor','interp','EdgeColor','interp');
%title(sprintf('Time = %.0f',tt{3}),'fontsize',12,'fontname','times')
%view(25,20);
% axis vis3d;
%axis([0 LX(3) 0 LZ(3) -6 -2]);

% hold on;
% slice(xx,zz,yy,permute(omega_x,[3,1,2]),[],LZ(3),[]);shading interp
% slice(xx,zz,yy,permute(omega_x,[3,1,2]),0,[],[]);shading interp
% slice(xx,zz,yy,permute(omega_x,[3,1,2]),[],[],-5);shading interp

% slice(xx,zz,yy,permute(omega_x,[3,1,2]),[],[],-4);shading flat

% colormap(flipud(cbrewer('div', 'RdYlBu', 100)));
% 
% % colormap(flipud(cm_1));
% xlabel('X','fontsize',12);
% ylabel('Y','fontsize',12);
% zlabel('Z','fontsize',12);
% alpha(1);
% box on; grid on
% caxis([-6 6]*1e-1);
% % co=colorbar('position',[.37 0.55 .015 .25]);
% co=colorbar;
% ylabel(co,'\omega_{x}','fontsize',12);
% 
% print -djpeg 3D_omega_16.jpg




%% plot the isosurfaces in terms of TH1 with different colors
% clear p layer color
% figure('position',[50,50,800,600]);
% i=1;
% k=2;
% layer = .6:.1:1.4;
% % color = jet(length(layer));
% cm=cmocean('thermal');
% color=interp1(1:length(cm),cm,linspace(1,length(cm),length(layer)));
% for j=1:length(layer)
% p(j) = patch(isosurface(xx(1:k:end,1:k:end,1:k:end),zz(1:k:end,1:k:end,1:k:end),yy(1:k:end,1:k:end,1:k:end),permute(TH1(1:k:end,1:k:end,1:k:end,i),[3,1,2]),layer(j)));
% isonormals(xx(1:k:end,1:k:end,1:k:end),zz(1:k:end,1:k:end,1:k:end),yy(1:k:end,1:k:end,1:k:end),permute(TH1(1:k:end,1:k:end,1:k:end,i),[3,1,2]),p(j));
% p(j).FaceColor = color(j,:);
% p(j).EdgeColor = 'none';
% hold on;
% end
% 
% view(25,20);
% camlight 
% % lighting gouraud
% axis([0 LX(3) 0 LZ(3) -5 -2]);
% alpha(1);
% box on; 
% 
% 
% xlabel('X','fontsize',12);
% ylabel('Y','fontsize',12);
% zlabel('Z','fontsize',12);
% alpha(1);
% box on; grid on

%%%
% print -djpeg 3D_TH1_16.jpg

%% multiple slices
% figure;
% % slice(xx,zz,yy,permute(omega_x,[3,1,2]),0:5:15.7,[],[]);shading flat
% % hold on;
% contourslice(xx,zz,yy,permute(TH1,[3,1,2]),[],1:40:120,[],50);%shading flat
% xlabel('X','fontsize',12);
% ylabel('Y','fontsize',12);
% zlabel('Z','fontsize',12);
% axis([0 LX 0 LZ -5 0]);
% alpha(.9);
% box on; 
% title('$Re$=1000, $Pr$=1, $Ri$=0.12, (384x192x246), $LY$=120, $d$=3','fontname','times','fontsize',12,...
%     'fontweight','normal','interpreter','latex');
% print -djpeg Slices_LY120_d3_180.jpg

%% plot the averaged omega_x^2 and K3D in x-z plane
% omega_x_sq_mean = mean(omega_x.^2,3);
% close all;
% % filename=[base_dir '/movie.h5'];
% % plot_xy = VideoWriter('omega_K3D_16_1_3_0.5');
% % plot_xy.FrameRate = 20; plot_xy.Quality = 100;
% % open(plot_xy);
% 
% figure('position',[50,50,500,700]);
% for k=1%:length(timestep)
% ax1 = axes('position',[0.1 0.55 .7 .4]);
% 
% % [cc,hh] = contourf(x,y,omega_x_mean');
% [cc,hh] = contourf(x,y,squeeze(omega_x_sq_mean(:,:,k))',100);
% set(hh,'edgecolor','none');
% caxis([0 max(max(omega_x_sq_mean))*.5]);
% colormap(ax1,cbrewer('div', 'PuOr', 100));
% colorbar('position',[0.85,0.55 0.015 0.4]);
% text(1.5,y(1)*.8,'$\langle \omega_{x}^{2} \rangle_{y}$','interpreter','latex',...
% 'fontsize',18,'color','w');
% % xlabel('X'); 
% ylabel('Z');
% set(gca,'fontsize',12,'xticklabel',[],'tickdir','out',...
%     'ticklength',[0.015 0.015],'xminortick','on','yminortick','on');
% box on;
% title('$Re$=1000, $Pr$=1, $Ri$=0.16, $d$=3','fontname','times','fontsize',12,...
%     'fontweight','normal','interpreter','latex');
% ylim([y(1) 0])
% ax2 = axes('position',[0.1 0.12 .7 .4]);
% % pcolor(x,y,tke_3D_m');shading flat;axis equal tight
% [cc,hh] = contourf(x,y,squeeze(K3D{i}(:,:,k))',100);set(hh,'edgecolor','none');
% % hold on; [cc,hh]=contour(x,y,tke_3D_m',10,'color','k');%set(cc,'color','k');
% colormap(ax2,flipud(cbrewer('div', 'RdYlBu', 100)));
% caxis([0 max(max(K3D{i}))*.7]);
% colorbar('position',[0.85,0.12 0.015 0.4]);
% text(1.5,y(1)*.8,'$K_{3D}$','interpreter','latex',...
% 'fontsize',18,'color','w');
% title(sprintf('t=%.0f',tt{i}),'fontweight','normal','fontname','times','fontsize',8);
% xlabel('X'); ylabel('Z');
% set(gca,'fontsize',12,'tickdir','out','ticklength',[0.015 0.015],...
%     'xminortick','on','yminortick','on');box on;
% ylim([y(1) 0]);
% drawnow
% 
% % M=getframe(gcf);
% %  writeVideo(plot_xy,M);
% 
%  
% % clf;
% %  close
%  end

% close(plot_xy)

% title('$Re$=1000, $Pr$=1, $Ri$=0.12, (384x96x342), $LY$=30, $d$=7','fontname','times','fontsize',12,...
%     'fontweight','normal','interpreter','latex');
% print -djpeg omegax_K3D_0.16_d3_252.jpg

%% PSD to estimate the wavelength explicitly 
% d=3;        % the distance to the boundary
% po = find(abs(y-(d-LY/2))==min(abs(y-(d-LY/2))));   % find the depth of shear layer
% dz = mean(diff(z)); kN = 2*pi/(2*dz);         % Nyquist wavenumber
% Sp_s=0;clear yn
% Coh_h = 0; Coh_l = 0;
% for j=1:size(U,2)
% for i=1:size(U,1)
%     yn(i,:) = (squeeze(omega_x(i,j,:)));              % the data
%     [Sp,kj]=fft_psd(yn(i,:)',dz,'rec');
%     Sp_s = Sp_s+Sp;                          % add all PSD together
% % confidence interval 
% DOF=2;         % no band-average to smooth the PSD
% Pro = .95; alpha = 1-Pro;
% q05=chi2inv(alpha/2,DOF); qF=chi2inv(1-alpha/2,DOF);
% Coh_h = Coh_h+DOF/(q05)*Sp; Coh_l = Coh_l+DOF/(qF)*Sp;    % this is in MATLAB definition
% 
% end
% end
% Coh_h = Coh_h/i/j; Coh_l = Coh_l/i/j; 
% Sp = Sp_s./i/j;                              % take average 
% [~,popo] = sort(Sp,'descend');
% 
% 
% % calculate the maximum wavenumber in radian
% L_max1 = 1./(kj(popo(1)));
% L_max2 = 1./(kj(popo(2)));
% L_max4 = 1./(kj(popo(4)));
% kj_max = (kj(popo(1))).*2*pi;
% kj_max2 = (kj(popo(2))).*2*pi;
% kj_max4 = (kj(popo(4))).*2*pi;
% 
%  
% 
% close all;
% figure('position',[50 50 700 800]);
% axes('position',[0.15 0.6 0.8 0.35]);
% hold on;plot(z,mean(squeeze(mean(omega_x,1))),'-k','linewidth',1.5); xlim([0 z(end)]);
% % plot(z,squeeze(omega_x(size(U,1)/2-10,po(1)+10,:)),'-k','linewidth',1.5); xlim([0 z(end)]);
% plot([0,z(end)],[0 0],'-k','linewidth',1.2);
% xlabel('Y','fontsize',12,'fontname','times'); 
% ylabel({'$\langle{\omega_{x}}\rangle$'},'fontsize',12,'fontname','times',...
%     'interpreter','latex'); 
% 
% set(gca,'fontsize',12);
% title('$Re$=1000, $Pr$=1, $Ri$=0.12, (345x192x343), $LY$=120, $d$=3','fontname','times','fontsize',12,...
%     'fontweight','normal','interpreter','latex');
% % title('$Re$=1000, $Pr$=1, $Ri$=0.12, (196x375x125), $LY$=30, $d$=3','fontname','times','fontsize',12,...
% %     'fontweight','normal','interpreter','latex');
% 
% 
% box on;
% axes('position',[0.15 0.12 0.8 0.4]);
% loglog(2*pi*kj,Sp,'k-','linewidth',1.1);hold on;
% % loglog([L_max2,L_max2],[1e-16,1e0],'k--')
% loglog([kj_max,kj_max],[1e-16,1e0],'k--');
% text(kj_max*1.1,1e-3,sprintf('%s %.4f','\lambda_y=',L_max1),'fontsize',12);
% % loglog([kj_max4,kj_max4],[1e-16,1e0],'k--');
% % text(kj_max4*1.1,10^(-3.5),sprintf('%s %.4f','\lambda=',L_max4),'fontsize',12);
% hold on; 
% pos=30;
% errorbar(2*pi*kj(pos),Sp(pos)*10.^-.5,10.^-.5*(Sp(pos)-Coh_l(pos)),...
%     10.^-.5*(Coh_h(pos)-Sp(pos)),0,0,'.k','linewidth',.9,'markersize',15);
% text(2*pi*kj(pos)*1.02,10.^-.5*Sp(pos),'95%','fontsize',12);
% 
% % plot the -5/3 slope (turbulence inertial range)
% hold on; plot(15*kj,8e-4*kj.^(-5/3),'k--','linewidth',1.1);
% text(2e-1,.5,'$l^{-5/3}$','interpreter','latex','fontsize',12);
% 
% xlabel('Spanwise Wavenumber, $l$','interpreter','latex','fontsize',12);
% ylabel('$PSD$','interpreter','latex','fontsize',12);
% axis([0 kN 1e-4 1e0]);
% set(gca,'fontsize',12,'tickdir','in','ticklength',[0.015 0.015],...
%     'ytick',10.^[-16:1],'xtick',[10^-1,10^0,10^1]);
% text(6e-2,1e-3,'d.o.f. = 2','fontsize',12);
% 
% set(gca,'xminortick','on','yminortick','on')
% % 
% print -djpeg LY_120_d3_180.jpg
%% plot PSD with multiple cases together to compare
% to make the plot work, you have to load different cases separatively

% d=7;        % the distance to the boundary
% % po = find(abs(y-(d-LY/2))==min(abs(y-(d-LY/2))));   % find the depth of shear layer
% dz = mean(diff(z)); kN = 1/(2*dz);         % Nyquist wavenumber
% Sp_s=0;clear yn Sp
% 
% for j=1:size(U,2)
% for i=1:size(U,1)
%     yn(i,:) = (squeeze(omega_x(i,j,:)));              % the data
%     [Sp,kj]=fft_psd(yn(i,:)',dz,'rec');
%     Sp_s = Sp_s+Sp;                          % add all PSD together
% end
% end
% Sp = Sp_s./i/j;                              % take average 
% [~,popo] = sort(Sp,'descend');
% % calculate the wavenumber of the leading modes in radian
% L_max1 = 1./(kj(popo(1)));
% L_max2 = 1./(kj(popo(2)));
% L_max4 = 1./(kj(popo(4)));
% kj_max = (kj(popo(1))).*2*pi;
% kj_max2 = (kj(popo(2))).*2*pi;
% kj_max4 = (kj(popo(4))).*2*pi;
% 
% % figure('position',[50 50 600 400]);
% 
% hold on;
% loglog(2*pi*kj,Sp,'-r','linewidth',1.1);hold on;
% 
% xlabel('Spanwise Wavenumber','interpreter','latex','fontsize',12);
% ylabel('$PSD$','interpreter','latex','fontsize',12);
% % axis([0 2*pi*kN 1e-5 1e0]);
% axis([0 6 1e-4 1e0]);
% set(gca,'fontsize',12,'tickdir','in','ticklength',[0.02 0.02],...
%     'ytick',10.^[-16:1],'xtick',[10^-1,10^0,10^1]);
% set(gca,'xminortick','on','yminortick','on')
% title('$Re$=1000, $Pr$=1, $Ri$=0.12, $d$=7','fontsize',12,...
%     'fontweight','normal','interpreter','latex');


%% after subsequently running all cases at previous section
% le = legend('$LY$ = 80','$LY$ = 120','$LY$ = 160');
% set(le,'interpreter','latex');
% 
% % loglog([L_max2,L_max2],[1e-16,1e0],'k--')
% hold on;
% loglog([kj_max,kj_max],[1e-4,1e0],'k--');
% text(kj_max*1.1,1e-3,sprintf('%s %.4f','\lambda_y=',L_max1),'fontsize',12);
% % loglog([kj_max4,kj_max4],[1e-16,1e0],'k--');
% % text(kj_max4*1.1,10^(-3.5),sprintf('%s %.4f','\lambda=',L_max4),'fontsize',12);
% 
% print -djpeg LY_80_120_160_sym.jpg




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


%%
% axes('position',[0.5 0.55 .35 .35]);
% i= 2;
% [faces,verts,colors] = isosurface(xx,zz,yy,TH1(:,:,:,i),1,omega_x(:,:,:,i));
% patch('Vertices',verts,'Faces',faces,'FaceVertexCData',colors,...
%     'FaceColor','interp','EdgeColor','interp');
% title(sprintf('Time = %.0f',tt(i)),'fontsize',12,'fontname','times')
% view(25,43)
% axis vis3d;
% colormap(flipud(cm_1))
% xlabel('X','fontsize',12);
% ylabel('Y','fontsize',12);
% zlabel('Z','fontsize',12);
% axis([0 LX 0 LZ -4 0]);
% alpha(.5)
% box on; 
% caxis([-2 2]*1e-3);
% co=colorbar('position',[.85 0.55 .015 .25]);
% ylabel(co,'\omega_{x}','fontsize',12);
% 
% axes('position',[0.02 0.1 .35 .35]);
% i= 3;
% [faces,verts,colors] = isosurface(xx,zz,yy,TH1(:,:,:,i),1,omega_x(:,:,:,i));
% patch('Vertices',verts,'Faces',faces,'FaceVertexCData',colors,...
%     'FaceColor','interp','EdgeColor','interp');
% title(sprintf('Time = %.0f',tt(i)),'fontsize',12,'fontname','times')
% view(25,43)
% axis vis3d;
% colormap(flipud(cm_1))
% xlabel('X','fontsize',12);
% ylabel('Y','fontsize',12);
% zlabel('Z','fontsize',12);
% axis([0 LX 0 LZ -4 0]);
% alpha(.5)
% box on; 
% caxis([-2 2]*1e-2);
% co=colorbar('position',[.37 0.1 .015 .25]);
% ylabel(co,'\omega_{x}','fontsize',12);
% 
% axes('position',[0.5 0.1 .35 .35]);
% i= 4;
% [faces,verts,colors] = isosurface(xx,zz,yy,TH1(:,:,:,i),1,omega_x(:,:,:,i));
% patch('Vertices',verts,'Faces',faces,'FaceVertexCData',colors,...
%     'FaceColor','interp','EdgeColor','interp');
% title(sprintf('Time = %.0f',tt(i)),'fontsize',12,'fontname','times')
% view(25,43)
% axis vis3d;
% colormap((cm_1));
% % cmocean('deep')
% xlabel('X','fontsize',12);
% ylabel('Y','fontsize',12);
% zlabel('Z','fontsize',12);
% axis([0 LX 0 LZ -4 0]);
% alpha(.8)
% box on; 
% co=colorbar('position',[.85 0.1 .015 .25]);
% ylabel(co,'\omega_{x}','fontsize',12);
% caxis([-2 2]*1e-2);
