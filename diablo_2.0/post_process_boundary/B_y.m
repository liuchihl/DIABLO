% this script is for calculating Rayleigh number, in which we need the
% spanwise averaged buoyancy field


clear;%clc; 
addpath(genpath('/glade/u/home/liuchihl/TOOLS'));
% input coefficients
for cc=[1:6]
if cc==1
%1) d=10, Ri=0.12
  fname1 = 'd_10/boundary/Ri_0.12_0.1_'; 
  RI=0.12; Re=1000; Pr=1;
  LX=28.28; NX=512;
  LY=20;    NY=361; 	 
  LZ=7.07;  NZ=128;
  load('TKEspectra_subKH_12_10_sta.mat','TKE','Time');
  num=10;
elseif cc==2
%2) d=6, Ri=0.12
  fname1 = 'd_6/boundary/Ri_0.12_0.1_';
  RI=0.12; Re=1000; Pr=1;
  LX=28.56; NX=512;
  LY=20;    NY=361;
  LZ=7.14;  NZ=128;
  load('TKEspectra_subKH_12_6_sta.mat','TKE','Time');
  num=10;
elseif cc==3
%3) d=4, Ri=0.12
  fname1 = 'd_4/boundary/Ri_0.12_0.1_';
  RI=0.12; Re=1000; Pr=1;
  LX=28.21; NX=512;
  LY=20;    NY=361;
  LZ=7.05;  NZ=128;
  load('TKEspectra_subKH_12_4_sta.mat','TKE','Time');
  num=15;
elseif cc==4
%3) d=3, Ri=0.12
  fname1 = 'd_3/boundary/Ri_0.12_0.1_';
  RI=0.12; Re=1000; Pr=1;
  LX=28.36; NX=512;
  LY=20;    NY=361;
  LZ=7.09;  NZ=128;
  load('TKEspectra_subKH_12_3_sta.mat','TKE','Time');
  num=10;
elseif cc==5
%4) d=2.5, Ri=0.12
  fname1 = 'd_2.5/Ri_0.12_0.1_';
  RI=0.12; Re=1000; Pr=1;
  LX=27.93; NX=512;
  LY=20;    NY=361;
  LZ=6.98;  NZ=128;
  load('TKEspectra_subKH_12_025_sta.mat','TKE','Time');
  num=6;
else
%5) d=2.5, Ri=0.12
  fname1 = 'd_2/Ri_0.12_0.1_';
  RI=0.12; Re=1000; Pr=1;
  LX=29.16; NX=512;
  LY=20;    NY=361;
  LZ=7.29;  NZ=128;
  load('TKEspectra_subKH_12_2_sta.mat','TKE','Time');
  num=6;
end
    % find t2d here
%1) find the first peak of Kkh first for all ensemble cases
% num=6;
tkh=zeros(1,num);
for i=1:num
%    if cc==1 && i==3
%    t2d(i) = 128.44;
%    elseif cc==1 && i==4
%    t2d(i) = 129.386;
%    elseif cc==1 && i==7
%    t2d(i) = 115;
%    elseif cc==3 && i==3
%    t2d(i) = 141;
%    elseif cc==3 && i==5
%    t2d(i) = 132.8;
%    elseif cc==3 && i==10
%    t2d(i) = 155;
%    elseif cc==4 && i==5
%    t2d(i) = 134;
%    elseif cc==4 && i==9
%    t2d(i) = 139;
%    elseif cc==5 && i==2
%    t2d(i) = 190;
%    else
    [pks,locs] = findpeaks(TKE{i}(2,:), 'MinPeakHeight',3e-5);
    t = locs(1);
    tkh(i) = Time{i}(t);

   
end
if cc==1
    tkh10 = tkh;
elseif cc==2
    tkh6 = tkh;
elseif cc==3
    tkh4 = tkh;
elseif cc==4
    tkh3 = tkh;
elseif cc==5
    tkh25 = tkh;
else cc==6
    tkh2 = tkh;
end

%2) load the time of the 3d slices
base_dir=['/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/'];
timestep=[4000:400:10000];
for i=1:num
   tt{i}=zeros(size(timestep));

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
    filename=[base_dir, fname1, num2str(i) '/output/'];
    info = h5info([filename,timename]);
    %time
    tt{i}(k) = info.Groups.Attributes.Value;
    end
end

%3) find the timestep that is closest to tkh
for i=1:num
    [~,ind] = min(abs(tt{i}-tkh(i)));
    ts(i) = timestep(ind);      % the timestep that is close to tkh
end

for i=1:num
    i
B_z{i}=zeros(NX,NY);
% tt{m}=zeros(size(timestep));
base_dir=['/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/'];

filename=[base_dir, fname1, num2str(i) '/output/'];
filename_mean=[base_dir , fname1, num2str(i) '/mean.h5'];

  if (ts(i)<10)
    timename=['out0000' int2str(ts(i)) '.h5'];
  elseif (ts(i)<100)
    timename=['out000' int2str(ts(i)) '.h5'];
  elseif (ts(i)<1000)
    timename=['out00' int2str(ts(i)) '.h5'];
  elseif (ts(i)<10000)
    timename=['out0' int2str(ts(i)) '.h5'];
  else
    timename=['out' int2str(ts(i)) '.h5'];
  end

varname=['/Timestep/TH1'];

if exist([filename,timename])
   B = h5read([filename,timename],varname);
   B_z{i} = mean(B,3);
% info = h5info([filename,timename]);
%time
% tt{m}(k) = info.Groups.Attributes.Value;

varname=['/gyf/' '0001'];             % Y-COORDINATE
gyf(:)=h5read(filename_mean,varname);

% NX=size(U,1); NZ=size(U,3); NY=size(U,2);
x=linspace(0,LX,NX); y=gyf;%linspace(-LY(i)/2,LY(i)/2,NY); 
z=linspace(0,LZ,NZ);

end
end


if cc==1
    save B_yavg_10 B_z y x tt tkh10
elseif cc==2
    save B_yavg_6 B_z y x tt tkh6
elseif cc==3
    save B_yavg_4 B_z y x tt tkh4
elseif cc==4
    save B_yavg_3 B_z y x tt tkh3
elseif cc==5
    save B_yavg_25 B_z y x tt tkh25
else 
    save B_yavg_2 B_z y x tt tkh2
end

end
