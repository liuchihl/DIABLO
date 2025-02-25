% This script illustrates how to load in 3D model data and calculate 
% some basic diagnostics
% Run after readmean.m

clear;%clc; 
addpath(genpath('/glade/u/home/liuchihl/TOOLS'));
% input coefficients
for cc=[2,3]
if cc==1
%1) d=10, Ri=0.12
  fname1 = 'd_10/boundary/Ri_0.12_0.1_'; 
  RI=0.12; Re=1000; Pr=1;
  LX=28.28; NX=512;
  LY=20;    NY=361; 	 
  LZ=7.07;  NZ=128;
num=6;
elseif cc==2
%2) d=4, Ri=0.12
fname1 = 'd_4/boundary/Ri_0.12_0.1_';
RI=0.12; Re=1000; Pr=1;
LX=28.21; NX=512;
LY=20;    NY=361;
LZ=7.05;  NZ=128;
num=15;
elseif cc==3
%3) d=3, Ri=0.12
fname1 = 'd_3/boundary/Ri_0.12_0.1_';
RI=0.12; Re=1000; Pr=1;
LX=28.36; NX=512;
LY=20;    NY=361;
LZ=7.09;  NZ=128;
num=10;
elseif cc==4
%4) d=2.5, Ri=0.12
fname1 = 'd_2.5/Ri_0.12_0.1_';
RI=0.12; Re=1000; Pr=1;
LX=27.93; NX=512;
LY=20;    NY=361;
LZ=6.98;  NZ=128;
num=6;
else
%5) d=2, Ri=0.12
fname1 = 'd_2/Ri_0.12_0.1_';
RI=0.12; Re=1000; Pr=1;
LX=29.16; NX=512;
LY=20;    NY=361;
LZ=7.29;  NZ=128;
num=6;
end
%timestep=[1000:1000:20000];
% timestep=[21000:1000:60000];
timestep=[400:400:20000];
cm_1 = load('MPL_gnuplot.rgb');
 i = 1;
%num=6;
for m=1:num
    m
U=zeros(NX,NY,NZ);V=zeros(NX,NY,NZ);
W=zeros(NX,NY,NZ);TH1=zeros(NX,NY,NZ);

U3d=zeros(NX,NY,NZ);V3d=zeros(NX,NY,NZ);W3d=zeros(NX,NY,NZ);
TH3d=zeros(NX,NY,NZ); TH3d_rms{m}=zeros(NX,NY,length(timestep));
U2d=zeros(NX,NY);V2d=zeros(NX,NY);W2d=zeros(NX,NY);
K2d{m}=zeros(NX,NY,length(timestep));K3d{m}=zeros(NX,NY,length(timestep));
 U3d_x=zeros(NX,NY,NZ);U3d_y=zeros(NX,NY,NZ);U3d_z=zeros(NX,NY,NZ);
 V3d_y=zeros(NX,NY,NZ);V3d_x=zeros(NX,NY,NZ);V3d_z=zeros(NX,NY,NZ);
 W3d_x=zeros(NX,NY,NZ);W3d_y=zeros(NX,NY,NZ);W3d_z=zeros(NX,NY,NZ);
%TH3d_x=zeros(NX,NY,NZ); TH3d_y=zeros(NX,NY,NZ);
%TH3d_z=zeros(NX,NY,NZ);

 SP_background{m}=zeros(size(timestep));
 SP_2d{m}=zeros(size(timestep));
 BF3d{m}=zeros(size(timestep)); D3d{m}=zeros(size(timestep)); 
 SP_sheardeform{m}=zeros(size(timestep));

dTHzdy=zeros(NX,NY); dTHzdx=zeros(NX,NY);
LHS=zeros(size(timestep));
D3d_p=zeros(size(timestep));
tt{m}=zeros(size(timestep));
%if     i==1    
%     base_dir=['/glade/campaign/univ/uosc0024/diablo_2.0/Large_case/'];
base_dir=['/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/'];
%elseif i==2
%    base_dir=['/glade/campaign/univ/uosc0024/diablo_2.0/Large_case/',fname2,'output/'];
%else
%    %base_dir=['/glade/campaign/univ/uosc0024/diablo_2.0/Large_case/',fname3,'output/'];
%    base_dir=['/glade/campaign/univ/uosc0024/diablo_2.0/Large_case/',fname3,'output/'];
%
%end
%base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/KH_test4/output/';

filename=[base_dir, fname1, num2str(m) '/output/'];
filename_mean=[base_dir , fname1, num2str(m) '/mean.h5'];
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
tt{m}(k) = info.Groups.Attributes.Value;

varname=['/gyf/' '0001'];             % Y-COORDINATE
gyf(:)=h5read(filename_mean,varname);

% NX=size(U,1); NZ=size(U,3); NY=size(U,2);
x=linspace(0,LX,NX); y=gyf;%linspace(-LY(i)/2,LY(i)/2,NY); 
z=linspace(0,LZ,NZ);


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
K3d{m}(:,:,k) = mean(.5*(U3d.^2+V3d.^2+W3d.^2),3);  %
TH3d_rms{m}(:,:,k) = rms(TH3d,3);
%K3d_umez{i}(k) = mean(mean(mean(.5*((U-UME_Z).^2+().^2+W3d.^2))));
% 3D perturbation derivatives

 for j=1:NY
 
  U3d_x(:,j,:) = gradient(squeeze(U3d(:,j,:))',mean(diff(x)))';
  V3d_x(:,j,:) = gradient(squeeze(V3d(:,j,:))',mean(diff(x)))';
  W3d_x(:,j,:) = gradient(squeeze(W3d(:,j,:))',mean(diff(x)))';
 % TH3d_x(:,j,:)= gradient(squeeze(TH3d(:,j,:))',mean(diff(x)))';
 end
 
 for j=1:NX
 
  V3d_y(j,:,:) = gradient(squeeze(V3d(j,:,:))',mean(diff(y)))';
  U3d_y(j,:,:) = gradient(squeeze(U3d(j,:,:))',mean(diff(y)))';
  W3d_y(j,:,:) = gradient(squeeze(W3d(j,:,:))',mean(diff(y)))';
 %TH3d_y(j,:,:) = gradient(squeeze(TH3d(j,:,:))',mean(diff(y)))';
 
 
  U3d_z(j,:,:) = gradient(squeeze(U3d(j,:,:)),mean(diff(z)));
  V3d_z(j,:,:) = gradient(squeeze(V3d(j,:,:)),mean(diff(z)));
  W3d_z(j,:,:) = gradient(squeeze(W3d(j,:,:)),mean(diff(z)));
 % TH3d_z(j,:,:) = gradient(squeeze(TH3d(j,:,:)),mean(diff(z)));
 
 end

% calculate 2D perturbation kinetic energy (Smyth et al., 2005)
 U2d = mean(U,3)-mean(mean(U,1),3);
 V2d = mean(V,3)-mean(mean(V,1),3);
 W2d = mean(W,3)-mean(mean(W,1),3);
 K2d{m}(:,:,k) = mean(.5*(U2d.^2+V2d.^2+W2d.^2),3);  %

 % % calculate the K3d energy budget
  % shear production: extract energy from the mean background field
 
  SP_background{m}(k) = -mean(mean(mean(...
      U3d.*V3d.*repmat(gradient(mean(mean(U,1),3),mean(diff(y))),NX,1,NZ) ))); 
  		
  % shear production: extract energy from the 2D KH flow field
  SP_2d{m}(k) = -mean(mean(mean(...
      U3d.^2.*repmat(gradient(U2d',mean(diff(x)))',1,1,NZ)+...
      V3d.^2.*repmat(gradient(V2d,mean(diff(y))),1,1,NZ) )));
 
  % shearing deformation
  SP_sheardeform{m}(k) = -mean(mean(mean(...
      U3d.*V3d.*(repmat(gradient(U2d,mean(diff(y))),1,1,NZ)+repmat(gradient(V2d',mean(diff(x)))',1,1,NZ) ) )));
 		
  % buoyancy flux
  BF3d{m}(k) = RI*mean(mean(mean(...
      TH3d.*V3d)));
  % dissipation
  D3d{m}(k) = -1/Re*mean(mean(mean(...
       2*(U3d_x.^2+V3d_y.^2+W3d_z.^2)+U3d_y.^2+U3d_z.^2+...
       V3d_x.^2+V3d_z.^2+W3d_x.^2+W3d_y.^2+...
       2*(U3d_z.*W3d_x+U3d_y.*V3d_x+V3d_z.*W3d_y) )));
k

% calculate the buoyancy variance budget equation
% dissipation rate (epsilon_p) 
%dTHzdy{i} = gradient(squeeze(mean(TH1,3)),mean(diff(y)));  
%dTHzdx{i} = gradient(squeeze(mean(TH1,3))',mean(diff(x)))';   

%D3d_p(k) = -1/(Re*Pr)*mean(mean(mean(...
%     (TH3d_x.^2+TH3d_y.^2+TH3d_z.^2) )));
% buoyancy flux
% BF3d_v{i}(k) = -mean(mean(mean(...
%     dTHzdy{i}.*TH3d.*V3d )));
% BF3d_h{i}(k) = -mean(mean(mean(...
%     dTHzdx{i}.*TH3d.*U3d )));
%LHS{i}(k) = mean(mean(mean(...
%    1/2*TH3d.^2 )));
%end
end
end
if cc==1
save K3D_budget_12_10 SP_background SP_2d SP_sheardeform BF3d D3d tt
 save ThreeD_xz12_10 K2d K3d TH3d_rms tt
elseif cc==2
save K3D_budget_12_4 SP_background SP_2d SP_sheardeform BF3d D3d tt
 save ThreeD_xz12_4 K2d K3d TH3d_rms tt
elseif cc==3
 save K3D_budget_12_3 SP_background SP_2d SP_sheardeform BF3d D3d tt
 save ThreeD_xz12_3 K2d K3d TH3d_rms tt
elseif cc==4
 save K3D_budget_12_25 SP_background SP_2d SP_sheardeform BF3d D3d tt
 save ThreeD_xz12_25 K2d K3d TH3d_rms tt
else 
 save K3D_budget_12_2 SP_background SP_2d SP_sheardeform BF3d D3d tt
 save ThreeD_xz12_2 K2d K3d TH3d_rms tt
end

end
%save B3d_budget_d10 D3d_p BF3d_v BF3d_h LHS tt

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

