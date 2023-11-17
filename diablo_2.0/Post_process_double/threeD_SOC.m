clear;
% D=2, Ri=0.16
fname1 = 'doubleshearlayer/Ri016/R1/3D/D_2_';
RI=0.16; Re=1000; Pr=1;
LX=29.92; NX=576;
LY=30;    NY=613;
LZ=7.48;  NZ=144;
m=1;

fname1 = 'doubleshearlayer/Ri016/R1/3D/D_3_';
 RI=0.16; Re=1000; Pr=1;
 LX=28.56; NX=576;
 LY=30;    NY=613;
 LZ=7.14;  NZ=144;
 m=2;
timestep = [50,55,60,65].*400;
%timestep = 50*400;
base_dir=['/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/'];
% double shear layer
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

if exist([filename,timename])
Ua(:,:,:)=h5read([filename,timename],varname1);
Va(:,:,:)=h5read([filename,timename],varname2);
Wa(:,:,:)=h5read([filename,timename],varname3);
%TH1(:,:,:)=h5read([filename,timename],varname4);

p2 = 409; p1 = 205;   % z = -5 ~ 5
%p2 = 358; p1 = 256;   % z = -2.5 ~ 2.5
U = Ua(:,p1:p2,:); V = Va(:,p1:p2,:); W = Wa(:,p1:p2,:);

info = h5info([filename,timename]);
%time
tt(k) = info.Groups.Attributes.Value;

varname=['/gyf/' '0001'];             % Y-COORDINATE
gyf(:)=h5read(filename_mean,varname);

% NX=size(U,1); NZ=size(U,3); NY=size(U,2);
x=linspace(0,LX,NX); y=gyf(p1:p2);
z=linspace(0,LZ,NZ);
end

% calculate tke
Up = U-mean(U,[1,3]);Vp = V-mean(V,[1,3]);Wp = W-mean(W,[1,3]);
tke(:,:,:,k) = .5.*(Up.^2+Vp.^2+Wp.^2);

% calculate epsilon

 for j=1:size(U,2)

  U_x(:,j,:) = gradient(squeeze(Up(:,j,:))',mean(diff(x)))';
  V_x(:,j,:) = gradient(squeeze(Vp(:,j,:))',mean(diff(x)))';
  W_x(:,j,:) = gradient(squeeze(Wp(:,j,:))',mean(diff(x)))';
 end

 for j=1:NX
  V_y(j,:,:) = gradient(squeeze(Vp(j,:,:))',mean(diff(y)))';
  U_y(j,:,:) = gradient(squeeze(Up(j,:,:))',mean(diff(y)))';
  W_y(j,:,:) = gradient(squeeze(Wp(j,:,:))',mean(diff(y)))';

  U_z(j,:,:) = gradient(squeeze(Up(j,:,:)),mean(diff(z)));
  V_z(j,:,:) = gradient(squeeze(Vp(j,:,:)),mean(diff(z)));
  W_z(j,:,:) = gradient(squeeze(Wp(j,:,:)),mean(diff(z)));
 end

  % dissipation
  eps(:,:,:,k) = -1/Re*2*((U_x.^2+V_y.^2+W_z.^2)+U_y.^2+U_z.^2+...
       V_x.^2+V_z.^2+W_x.^2+W_y.^2+...
       2*(U_z.*W_x+U_y.*V_x+V_z.*W_y));


% energy containing scale
L(:,:,:,k) = tke(:,:,:,k).^3./eps(:,:,:,k);

end

save L_en_D3.mat L eps tke tt

