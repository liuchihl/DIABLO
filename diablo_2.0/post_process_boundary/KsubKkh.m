% this script is to calculate the subharmonic TKE spectru,m Ksub andKH TKE spectrum, K_KH
clear;
addpath(genpath('/glade/u/home/liuchihl/TOOLS'));

% Input coefficients in DIABLO coordinate
for cc=[2]
	if cc==1
%1) 0.12_10
RI=0.12; Re=1000;Pr=1;
LX=28.28; NX=512;
LY=20;    NY=361;
LZ=7.07;  NZ=128;
base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_10/boundary/Ri_0.12_0.1_';
num=10;
	elseif cc==2

%2) 0.12_6
 RI=0.12; Re=1000;Pr=1;
 LX=28.56; NX=512;
 LY=20;    NY=361;
 LZ=7.14;  NZ=128;
 base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_6/boundary/Ri_0.12_0.1_';
num=10;
	elseif cc==3
%3) 0.12_4
 RI=0.12; Re=1000;Pr=1;
 LX=28.21; NX=512;
 LY=20;    NY=361;
 LZ=7.05;  NZ=128;
 base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_4/boundary/Ri_0.12_0.1_';
num=15;
%4) 0.12_3
	elseif cc==4
RI=0.12; Re=1000;Pr=1;
LX=28.36; NX=512;
LY=20;    NY=361;
LZ=7.09;  NZ=128;
base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_3/boundary/Ri_0.12_0.1_';
num=10;
%5) 0.12_2
	elseif cc==5
RI=0.12; Re=1000;Pr=1;
LX=29.16; NX=512;
LY=20;    NY=361;
LZ=7.29;  NZ=128;
base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_2/Ri_0.12_0.1_';
num=6;
%5) 0.16_1
%	elseif cc==5
%RI=0.16; Re=1000;Pr=1;
%LX=28.64; NX=512;
%LY=20;    NY=361;
%LZ=7.16;  NZ=128;
%base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_3/boundary/Ri_0.16_0.1_';
	else
RI=0.12; Re=1000;Pr=1;
LX=27.93; NX=512;
LY=20;    NY=361;
LZ=6.98;  NZ=128;
base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_2.5/Ri_0.12_0.1_';
num=6;
	end
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

if cc==3
	if i~=2
	k1=1:nk;
	else i==2
	k1=1002:1502;
	end
else 
	k1=1:nk;
end
for k=k1
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
%varname=['/time/' timename];            % TIME
%time(k)=h5read(filename_mean,varname);

varname=['/u_xz/' timename];
U=h5read(filename,varname);
Umean=mean(U,[1,2]); u=U-Umean;

varname=['/v_xz/' timename];
V=h5read(filename,varname);
Vmean=mean(V,[1,2]); v=V-Vmean;

varname=['/w_xz/' timename];
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
% Sp_ku(:,k) = Sp_su'./NZ;                          % take average
% Sp_kv(:,k) = Sp_sv'./NZ;
% Sp_kw(:,k) = Sp_sw'./NZ;

% if we only want the time series of the 2 lowest wavenumbers

Sp_ku(:,k) = Sp_su(1:2,:)./NZ;                          % take average
Sp_kv(:,k) = Sp_sv(1:2,:)./NZ;
Sp_kw(:,k) = Sp_sw(1:2,:)./NZ;

if  cc==3 && i==2
Sp_ku(:,k-1001) = Sp_su(1:2,:)./NZ;                          % take average
Sp_kv(:,k-1001) = Sp_sv(1:2,:)./NZ;
Sp_kw(:,k-1001) = Sp_sw(1:2,:)./NZ;
end

end
TKE{i} = .5*(Sp_ku+Sp_kv+Sp_kw);
Time{i} = time;
end
if cc==1
save TKEspectra_subKH_12_10_sta TKE Time
elseif cc==2
save TKEspectra_subKH_12_6_sta TKE Time
elseif cc==3
save TKEspectra_subKH_12_4_sta TKE Time
elseif cc==4
save TKEspectra_subKH_12_3_sta TKE Time
elseif cc==5
save TKEspectra_subKH_12_2_sta TKE Time
%elseif cc==6
%save TKEspectra_subKH_16_3_sta TKE Time
else cc==6
save TKEspectra_subKH_12_025_sta TKE Time
end


end 
% save TKEspecta_subKH_12_4 TKE time
% save TKEspecta_subKH_16_10 TKE time
% save TKEspecta_subKH_16_4 TKE time

