% this script is to calculate the PSD of turbulent kinetic energy
clear;
addpath(genpath('/glade/u/home/liuchihl/TOOLS'));

% Input coefficients in DIABLO coordinate
for cc=[1:2]
	if cc==1
%1) 0.16_
RI=0.16; Re=1000;Pr=9;
LX=36.96; NX=2048;
LY=30;    NY=1657;
LZ=1;     NZ=1;
base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/doubleshearlayer/Ri016/R3/D_3_Pr9_res1';
num=1;
	elseif cc==2

%2) 0.12_6
 RI=0.16; Re=1000;Pr=9;
 LX=36.96; NX=1296;
 LY=30;    NY=1081;
 LZ=1;  NZ=1;
base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/doubleshearlayer/Ri016/R3/D_3_Pr9_res3';
num=1;
	end
clear Sp_ku Sp_kv Sp_kb TKE k2d
for i=1:num
filename=[base_dir, '/movie.h5'];
filename_mean=[base_dir, '/mean.h5'];
file_info=h5info(filename_mean);
att_info=file_info.Groups.Attributes;
nk=att_info.Value;

x=linspace(0,LX,NX);
z=linspace(0,LZ,NZ);
U=zeros(NX,NZ);V=zeros(NX,NZ);W=zeros(NX,NZ);
u=zeros(NX,NZ);v=zeros(NX,NZ);w=zeros(NX,NZ);

time=zeros(1,nk); 
TKE{i} = zeros(2,nk);
%time = zeros(1,nk);
k2d = zeros(NY,nk);
%TKE = zeros(length(kj),2);
	for k=1:nk

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
varname=['/gyf/' timename];             % vertical grids
gyf=h5read(filename_mean,varname);
varname=['/k2d/' timename];
k2d(:,k)=h5read(filename_mean,varname);
%only for 3D
%varname=['/k3d/' timename];
%k3d(:,k)=h5read(filename_mean,varname);
	end
% find t2d and t3d
[~,inx2d] = max(mean(k2d)); 
t2d = time(inx2d);
%[~,inx3d] = max(mean(k3d));
%t3d = time(inx3d);
k1 = [inx2d,inx2d+100];
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

varname=['/u_xy/' timename];
U=h5read(filename,varname);
Umean=mean(U,[1]); u=U-Umean;

varname=['/v_xy/' timename];
V=h5read(filename,varname);
Vmean=mean(V,[1]); v=V-Vmean;

varname=['/th1_xy/' timename];
B=h5read(filename,varname);
Bmean=mean(B,[1]); b=B-Bmean;
% for 2D, comment W out
%varname=['/w_xz/' timename];
%W=h5read(filename,varname);
%Wmean=mean(W,[1,2]); w=W-Wmean;

% PSD 
dx = mean(diff(x)); kN = 2*pi/(2*dx);         % Nyquist wavenumber
dy = mean(diff(gyf));
Sp_su=0;Sp_sv=0;Sp_sw=0;Sp_sb=0;
for j=1:size(U,2)                             % loop over vertical direction 
     [Sp_u,kj]=fft_psd(u(:,j),dx,'rec');
     Sp_su = Sp_su+Sp_u;                      % add all PSD together
     [Sp_v,kj]=fft_psd(v(:,j),dy,'rec');
     Sp_sv = Sp_sv+Sp_v;
     [Sp_b,kj]=fft_psd(b(:,j),dy,'rec');
     Sp_sb = Sp_sb+Sp_b;
%     [Sp_w,kj]=fft_psd(w(:,j),dx,'rec');
%     Sp_sw = Sp_sw+Sp_w;
end
% if we want the time series of the entire spectrum
% Sp_ku(:,k) = Sp_su'./NZ;                          % take average
% Sp_kv(:,k) = Sp_sv'./NZ;
% Sp_kw(:,k) = Sp_sw'./NZ;

% if we want the spectra at t2d and t3d

Sp_ku(:,k) = Sp_su./NY;                          % take average
Sp_kv(:,k) = Sp_sv./NY;
Sp_kb(:,k) = Sp_sb./NY;
%Sp_kw(:,k) = Sp_sw./NZ;

end
TKE = .5*(Sp_ku+Sp_kv);%+Sp_kw);
ks = 2*pi*kj;
%EPS = 2/Re*ks'.^2.*TKE;
end

if cc==1
save PSD_tke_Pr9_Ri16_D_3_res1 TKE Sp_kb ks t2d k2d 
else cc==2
save PSD_tke_Pr9_Ri16_D_3_res3 TKE Sp_kb ks t2d k2d
end


end
