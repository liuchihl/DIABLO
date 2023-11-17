% This script computes PSD at t2d
% some basic diagnostics
% Run after readmean.m

clear;%clc;
addpath(genpath('/glade/u/home/liuchihl/TOOLS'));
% input coefficients
D = [0.5,1,2,3];
for cc=[2:4]

if cc==1
%1) D=0.5, Ri=0.16
  fname1 = 'doubleshearlayer/Ri016/R1/3D/D_0.5_';
  RI=0.16; Re=1000; Pr=1;
  LX=36.96; NX=768;
  LY=30;    NY=613;
  LZ=9.24;  NZ=192;
timestep = 14*400;
elseif cc==2
%2) D=1, Ri=0.16
fname1 = 'doubleshearlayer/Ri016/R1/3D/D_1_';
RI=0.16; Re=1000; Pr=1;
LX=78.54; NX=1536;
LY=30;    NY=613;
LZ=19.64;  NZ=384;
timestep = 44*400;
elseif cc==3
%3) D=2, Ri=0.16
fname1 = 'doubleshearlayer/Ri016/R1/3D/D_2_';
RI=0.16; Re=1000; Pr=1;
LX=29.92; NX=576;
LY=30;    NY=613;
LZ=7.48;  NZ=144;
timestep = 36*400;
elseif cc==4
%4) D=3, Ri=0.16
fname1 = 'doubleshearlayer/Ri016/R1/3D/D_3_';
RI=0.16; Re=1000; Pr=1;
LX=28.56; NX=576;
LY=30;    NY=613;
LZ=7.14;  NZ=144;
timestep = 38*400;
else cc==5
%5) D=infinity, Ri=0.16, butterfly case #2
fname1 = 'butterfly/Ri_0.16_0.05_';
RI=0.16; Re=1000; Pr=1;
LX=27.76; NX=512;
LY=20;    NY=361;
LZ=6.94;  NZ=128;
num=4;
end

base_dir=['/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/',fname1];
i=1;
filename=[base_dir, num2str(i) '/output/'];
filename_mean=[base_dir, num2str(i) '/mean.h5'];
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



U(:,:,:)=h5read([filename,timename],varname1);
V(:,:,:)=h5read([filename,timename],varname2);
W(:,:,:)=h5read([filename,timename],varname3);
TH1(:,:,:)=h5read([filename,timename],varname4);

info = h5info([filename,timename]);
%time
tt{i}(k) = info.Groups.Attributes.Value;

varname=['/gyf/' '0001'];             % Y-COORDINATE
gyf(:)=h5read(filename_mean,varname);

x=linspace(0,LX,NX); y=gyf; z=linspace(0,LZ,NZ);

Umean=mean(U,[1,3]); u=U-Umean;

Vmean=mean(V,[1,3]); v=V-Vmean;

Wmean=mean(W,[1,3]); w=W-Wmean;

% ind = find( abs(y-D(cc))==min(abs(y-D(cc)) ) );
% PSD
dx = mean(diff(x)); kN = 2*pi/(2*dx);         % Nyquist wavenumber
Sp_su=0;Sp_sv=0;Sp_sw=0;

for m=1:size(U,3)
% j=ind;
for j=1:size(U,2)                           % loop over y,z direction
     [Sp_u,kj]=fft_psd(u(:,j,m),dx,'rec',1024);
     Sp_su = Sp_su+Sp_u;                          % add all PSD together
     [Sp_v,kj]=fft_psd(v(:,j,m),dx,'rec',1024);
     Sp_sv = Sp_sv+Sp_v;
     [Sp_w,kj]=fft_psd(w(:,j,m),dx,'rec',1024);
     Sp_sw = Sp_sw+Sp_w;
end
end
% if we want the entire spectrum at t2d for each cases of D
 Sp_ku(:,k) = Sp_su'./(NZ*NY);                          % take average
 Sp_kv(:,k) = Sp_sv'./(NZ*NY);
 Sp_kw(:,k) = Sp_sw'./(NZ*NY);
end

TKE{i} = .5*(Sp_ku+Sp_kv+Sp_kw);
%Time{i} = time;
kj = kj*2*pi;

if cc==1
save PSD_t2d_16_0.5_3D.mat TKE tt kj
elseif cc==2
save PSD_t2d_16_1_3D.mat TKE tt kj
elseif cc==3
save PSD_t2d_16_2_3D.mat TKE tt kj
elseif cc==4
save PSD_t2d_16_3_3D.mat TKE tt kj
else cc==5
save PSD_t_16_butterfly TKE Time
end

clear;

end
