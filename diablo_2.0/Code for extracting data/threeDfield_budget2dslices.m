% This script illustrates how to load in 3D model data and calculate 
% some basic diagnostics 

clear;%clc; 
addpath(genpath('/glade/u/home/liuchihl/TOOLS'));
% input coefficients
 cc=1
if cc==1
%1) D=0.5, Ri=0.16
  fname1 = 'doubleshearlayer/Ri016/R1/3D/D_0.5_';
  RI=0.16; Re=1000; Pr=1;
  LX=36.96; NX=768;
  LY=30;    NY=613;
  LZ=9.24;  NZ=192;
num=2;
%timestep = 400*10;
  timestep = 400*[13,14,17];
%timestep = 400*[2:24];

elseif cc==2
%2) D=1, Ri=0.16
fname1 = 'doubleshearlayer/Ri016/R1/3D/D_1_';
RI=0.16; Re=1000; Pr=1;
LX=78.54; NX=1536;
LY=30;    NY=613;
LZ=19.64;  NZ=384;
num=1;
%timestep = 400*[49];  % max growth rate
timestep = 400*[44];  % earlier growth showing CCI
%timestep = 400*[30:50];
elseif cc==3
%3) D=2, Ri=0.16
fname1 = 'doubleshearlayer/Ri016/R1/3D/D_2_';
RI=0.16; Re=1000; Pr=1;
LX=29.92; NX=576;
LY=30;    NY=613;
LZ=7.48;  NZ=144; 
num=1;
timestep = 400*36;
elseif cc==4
%4) D=3, Ri=0.16
fname1 = 'doubleshearlayer/Ri016/R1/3D/D_3_';
RI=0.16; Re=1000; Pr=1;
LX=28.56; NX=576;
LY=30;    NY=613;
LZ=7.14;  NZ=144;
num=1;
%timestep = 400*30;
 timestep = 400*42;
% timestep = 400*[10:54];
else cc==5
%5) D=infinity, Ri=0.16, butterfly case #2
fname1 = 'butterfly/Ri_0.16_0.05_';
RI=0.16; Re=1000; Pr=1;
LX=27.76; NX=512;
LY=20;    NY=361;
LZ=6.94;  NZ=128;
num=4;
end

%X = linspace(0,LX,NX); Z = linspace(-LY/2,LY/2,NY);
%timestep=[4800,5600];  % d=10, case=3, at max KH amplitude

%timestep=[14800];
cm_1 = load('MPL_gnuplot.rgb');
% i = 1;
% num=1;
for m=num
    m
U=zeros(NX,NY,NZ);V=zeros(NX,NY,NZ);
W=zeros(NX,NY,NZ);TH1=zeros(NX,NY,NZ);

% U3d=zeros(NX,NY,NZ);V3d=zeros(NX,NY,NZ);W3d=zeros(NX,NY,NZ);
% TH3d=zeros(NX,NY,NZ); TH3d_rms{m}=zeros(NX,NY,length(timestep));
% U2d=zeros(NX,NY);V2d=zeros(NX,NY);W2d=zeros(NX,NY);
%K2d{m}=zeros(NX,NY,length(timestep));
% K3d=zeros(1,1);
% U3d_x=zeros(NX,NY,NZ);U3d_y=zeros(NX,NY,NZ);U3d_z=zeros(NX,NY,NZ);
% V3d_y=zeros(NX,NY,NZ);V3d_x=zeros(NX,NY,NZ);V3d_z=zeros(NX,NY,NZ);
% W3d_x=zeros(NX,NY,NZ);W3d_y=zeros(NX,NY,NZ);W3d_z=zeros(NX,NY,NZ);
% sh=zeros(NX,NY,length(timestep));
% st=zeros(NX,NY,length(timestep));
% uw=zeros(NX,NY,length(timestep));
% usws=zeros(NX,NY,length(timestep));
% bsh=zeros(NX,NY,length(timestep));
%TH3d_x=zeros(NX,NY,NZ); TH3d_y=zeros(NX,NY,NZ);
%TH3d_z=zeros(NX,NY,NZ);

%  SP_background=zeros(NX,NY,length(timestep));
%  SP_shear=zeros(NX,NY,length(timestep));
%  SP_strain=zeros(NX,NY,length(timestep));
%  BF3d=zeros(NX,NY,length(timestep));% D3d{m}=zeros(size(timestep)); 
% SP_sheardeform{m}=zeros(size(timestep));

%dTHzdy=zeros(NX,NY); dTHzdx=zeros(NX,NY);
%LHS=zeros(size(timestep));
%D3d_p=zeros(size(timestep));
omega_x = zeros(NX,NY,NZ,2);
tt=zeros(size(timestep));
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
tt(k) = info.Groups.Attributes.Value;

varname=['/gyf/' '0001'];             % Y-COORDINATE
gyf(:)=h5read(filename_mean,varname);

% NX=size(U,1); NZ=size(U,3); NY=size(U,2);
x=linspace(0,LX,NX); y=gyf;%linspace(-LY(i)/2,LY(i)/2,NY); 
z=linspace(0,LZ,NZ);

%U3d = U-mean(U,3);V3d = V-mean(V,3);W3d = W-mean(W,3);
% vorticities
% for j=1:size(U,3)
%wx(:,:,j) = mmderiv(x,squeeze(W(:,:,j))) ;
%vx(:,:,j) = mmderiv(x,squeeze(V(:,:,j))) ;
% wx3d(:,:,j) = mmderiv(x,squeeze(W3d(:,:,j,k))) ;
% vx3d(:,:,j) = mmderiv(x,squeeze(V3d(:,:,j,k))) ;
% end
% for j=1:size(U,3)
%wy(:,:,j) = mmderiv(y,squeeze(W(:,:,j))')' ;
%uy(:,:,j) = mmderiv(y,squeeze(U(:,:,j))')' ;
% wy3d(:,:,j) = mmderiv(y,squeeze(W3d(:,:,j,k))')' ;
% uy3d(:,:,j) = mmderiv(y,squeeze(U3d(:,:,j,k))')' ;
% end
% for j=1:size(U,2)
%uz(:,j,:) = mmderiv(z,squeeze(U(:,j,:))')' ;
%vz(:,j,:) = mmderiv(z,squeeze(V(:,j,:))')' ;
% uz3d(:,j,:) = mmderiv(z,squeeze(U3d(:,j,:,k))')' ;
% vz3d(:,j,:) = mmderiv(z,squeeze(V3d(:,j,:,k))')' ;
% end

%omega_x(:,:,:,k) = wy-vz;
%omega_y(:,:,:,k) = uz-wx;
%omega_z(:,:,:,k) = vx-uy;
% omega_x3d = wy3d-vz3d;
% omega_y3d = uz3d-wx3d;
% omega_z3d = vx3d-uy3d;
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
%U3dV3d = U3d.*V3d;
TH3d = TH1-mean(TH1,3);
%K3d(k) = mean(mean(mean(.5*(U3d.^2+V3d.^2+W3d.^2) )));  %
k3d(:,:,k) = mean( .5*(U3d.^2+V3d.^2+W3d.^2) ,3);  %
%TH3d_rms{m}(:,:,k) = rms(TH3d,3);
%K3d_umez{i}(k) = mean(mean(mean(.5*((U-UME_Z).^2+().^2+W3d.^2))));

% calculate 2D perturbation kinetic energy (Smyth et al., 2005)
 U2d = mean(U,3)-mean(mean(U,1),3);
 V2d = mean(V,3)-mean(mean(V,1),3);
 W2d = mean(W,3)-mean(mean(W,1),3);
  TH2d = mean(TH1,3)-mean(mean(TH1,1),3);
% K2d{m}(:,:,k) = mean(.5*(U2d.^2+V2d.^2+W2d.^2),3);  %

 % % calculate the K3d energy budget: but now we want the terms seperately

% 1) (du2d/dz+dw2d/dx): shear term
% du2dz(:,:,k) = gradient(U2d,mean(diff(y)));
% dw2dx(:,:,k) = gradient(V2d',mean(diff(x)))';
% sh(:,:,k) =  du2dz(:,:,k)+dw2dx(:,:,k);

% % 2)  (du2d/dx-dw2d/dz): strain term
% st(:,:,k) = mean( repmat(gradient(U2d',mean(diff(x)))',1,1,NZ)-repmat(gradient(V2d,mean(diff(y))),1,1,NZ) ,3);
% 
% % 3) (dU/dz+du2d/dz+dw2d/dx): background
% bsh(:,:,k) =  repmat(gradient(mean(mean(U,1),3),mean(diff(y))),NX,1) ; 
% 
% % 4) stresses: u3dw3d and -1/2(u3d^2-w3d^2)
%  uw(:,:,k) = mean(U3d.*V3d,3);
% usws(:,:,k) = -.5*mean(U3d.^2-V3d.^2,3);
% %uu(:,:,k) =mean(U3d.*U3d,3);
% %ww(:,:,k) =mean(V3d.*V3d,3);
% %vv(:,:,k) = mean(W3d.*W3d,3);
 % 5) shear production terms:
SP_shear(:,:,k) = ...
-mean(U3d.*V3d.*(repmat(gradient(U2d,mean(diff(y))),1,1,NZ)+repmat(gradient(V2d',mean(diff(x)))',1,1,NZ)),3);
%SP_shtime(k) = ...
%-mean(mean(mean(U3d.*V3d.*(repmat(gradient(U2d,mean(diff(y))),1,1,NZ)+repmat(gradient(V2d',mean(diff(x)))',1,1,NZ)) )));

SP_strain(:,:,k) =...
-.5*mean((U3d.^2-V3d.^2).*(repmat(gradient(U2d',mean(diff(x)))',1,1,NZ)-repmat(gradient(V2d,mean(diff(y))),1,1,NZ)),3);
SP_background(:,:,k) = ...
-mean(U3d.*V3d.*repmat(gradient(mean(mean(U,1),3),mean(diff(y))),NX,1,NZ),3);
% buoyancy flux
BF3d(:,:,k) = ...
RI*mean(TH3d.*V3d,3);
% bw(:,:,k) = mean(TH3d.*V3d,3);

%U_xy(:,k) = mean(mean(U,3),1);
% B_y(:,:,k) = mean(TH1,3);
%u2d(:,:,k) = U2d;
%w2d(:,:,k) = V2d;
%b2d(:,:,k) = TH2d;

%omega_x = mean(omega_x,3);
%omega_y = mean(omega_y,3);
%omega_z = mean(omega_z,3);


end
end
% save('K3dbudgetXZ_16_D3_1.mat','k3d','sh', 'st' ,'bsh', 'uw', 'usws','du2dz','dw2dx', 'SP_shear', 'SP_strain', 'SP_background', 'BF3d', 'U2d','V2d','tt', 'x' ,'y');
% save('sh_3_1.mat','sh','uw','tt', 'x' ,'y');
time = tt(1);
%save('vort_1_1_217.mat','omega_x','omega_y','omega_z','B_y','x','y','z','time');
%save('K3dbudgetXZ_16_D1_1_t_217.mat','k3d', 'SP_shear', 'SP_strain', 'SP_background', 'BF3d','tt', 'x' ,'y');
save('K3dbudgetXZ_16_D0.5_2_t.mat','k3d', 'SP_shear', 'SP_strain', 'SP_background', 'BF3d','tt', 'x' ,'y');

