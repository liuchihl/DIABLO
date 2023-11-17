% this script is to calculate the subharmonic TKE spectru,m Ksub andKH TKE spectrum, K_KH
clear;
addpath(genpath('/glade/u/home/liuchihl/TOOLS'));

% Input coefficients in DIABLO coordinate
for cc=[6]
if cc==1
%1) D=0.5, Ri=0.16
  fname1 = 'doubleshearlayer/Ri016/R1/3D/D_0.5_';
  RI=0.16; Re=1000; Pr=1;
  LX=36.96; NX=768;
  LY=30;    NY=613;
  LZ=9.24;  NZ=192;
num=1;
elseif cc==2
%2) D=1, Ri=0.16
fname1 = 'doubleshearlayer/Ri016/R1/3D/D_1_';
RI=0.16; Re=1000; Pr=1;
LX=78.54; NX=1536;
LY=30;    NY=613;
LZ=19.64;  NZ=384;
num=3;
elseif cc==3
%3) D=1.5, Ri=0.16
fname1 = 'doubleshearlayer/Ri016/R1/3D/D_1.5_';
RI=0.16; Re=1000; Pr=1;
LX=31.416; NX=576;
LY=30;    NY=613;
LZ=7.854;  NZ=144;
num=1;
elseif cc==4
%3) D=2, Ri=0.16
fname1 = 'doubleshearlayer/Ri016/R1/3D/D_2_';
RI=0.16; Re=1000; Pr=1;
LX=29.92; NX=576;
LY=30;    NY=613;
LZ=7.48;  NZ=144;
num=1;
elseif cc==5
%4) D=3, Ri=0.16
fname1 = 'doubleshearlayer/Ri016/R1/3D/D_3_';
RI=0.16; Re=1000; Pr=1;
LX=28.56; NX=576;
LY=30;    NY=613;
LZ=7.14;  NZ=144;
num=1;
else cc==6
%5) D=infinity, Ri=0.16, butterfly case #2
fname1 = 'butterfly/Ri_0.16_0.05_';
RI=0.16; Re=1000; Pr=1;
LX=27.76; NX=512;
LY=20;    NY=361;
LZ=6.94;  NZ=128;
num=4;
end
base_dir = ['/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/' fname1];
%num=15;
for i=1:num
filename=[base_dir, num2str(i) '/movie.h5'];
filename_mean=[base_dir num2str(i) '/mean.h5'];
file_info=h5info(filename_mean);
att_info=file_info.Groups.Attributes;
nk=att_info.Value;

x=linspace(0,LX,NX);
z=linspace(0,LZ,NZ);
U=zeros(NX,NZ);V=zeros(NX,NZ);W=zeros(NX,NZ);
u=zeros(NX,NZ);v=zeros(NX,NZ);w=zeros(NX,NZ);

time=zeros(1,nk); 
%Sp_ku = zeros(2,nk);Sp_kv = zeros(2,nk);
%Sp_kw = zeros(2,nk); TKE{i} = zeros(2,nk);
%time = zeros(1,nk);
%

for k=1:nk
k;
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
end

% eliminate time discontinuity result from multiple runs
k1 = 1:nk;
while any(diff(time)<0.1)
ind1 = find(diff(time)<0.1);
k1(ind1(1)+1) = [];
time(ind1(1)+1) = [];
end

for k=1:length(k1)

  if (k1(k)<10)
    timename=['000' int2str(k1(k))];
  elseif (k1(k)<100)
    timename=['00' int2str(k1(k))];
  elseif (k1(k)<1000)
    timename=['0' int2str(k1(k))];
  else
    timename=[int2str(k1(k))];
  end
%varname=['/time/' timename];            % TIME
%time(k)=h5read(filename_mean,varname);

%varname=['/u_xz/' timename]
varname=['/u_xy/' timename];;
U=h5read(filename,varname);
Umean=mean(U,[1,2]); u=U-Umean;

%varname=['/v_xz/' timename];;
varname=['/v_xy/' timename];
V=h5read(filename,varname);
Vmean=mean(V,[1,2]); v=V-Vmean;

%varname=['/w_xz/' timename];
varname=['/w_xy/' timename];
W=h5read(filename,varname);
Wmean=mean(W,[1,2]); w=W-Wmean;

% PSD 
dx = mean(diff(x)); kN = 2*pi/(2*dx);         % Nyquist wavenumber
Sp_su=0;Sp_sv=0;Sp_sw=0;
for j=1:size(U,2)                           % loop over z direction 
     [Sp_u,kj]=fft_psd(u(:,j),dx,'rec');
     Sp_su = Sp_su+Sp_u;                          % add all PSD together
     [Sp_v,kj]=fft_psd(v(:,j),dx,'rec');
     Sp_sv = Sp_sv+Sp_v;
     [Sp_w,kj]=fft_psd(w(:,j),dx,'rec');
     Sp_sw = Sp_sw+Sp_w;
end
% if we want the time series of the entire spectrum
 Sp_ku(:,k) = Sp_su'./NZ;                          % take average
 Sp_kv(:,k) = Sp_sv'./NZ;
 Sp_kw(:,k) = Sp_sw'./NZ;

% if we only want the time series of the 2 lowest wavenumbers

%Sp_ku(:,k) = Sp_su(1:2,:)./NZ;                          % take average
%Sp_kv(:,k) = Sp_sv(1:2,:)./NZ;
%Sp_kw(:,k) = Sp_sw(1:2,:)./NZ;


TKE{i} = .5*(Sp_ku+Sp_kv+Sp_kw);
Time{i} = time;
kj = kj*2*pi;
end
if cc==1
save PSD_t_16_0.5.mat TKE Time kj
elseif cc==2
save PSD_t_16_1.mat TKE Time kj
elseif cc==3
save PSD_t_16_1.5.mat TKE Time kj
elseif cc==4
save PSD_t_16_2.mat TKE Time kj
elseif cc==5
save PSD_t_16_3.mat TKE Time kj
else cc==6
save PSD_t_16_butterfly.mat TKE Time kj
end

end

%clear;
end 
% save TKEspecta_subKH_12_4 TKE time
% save TKEspecta_subKH_16_10 TKE time
% save TKEspecta_subKH_16_4 TKE time

