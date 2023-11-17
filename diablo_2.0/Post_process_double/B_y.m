% this script is for calculating Rayleigh number, in which we need the
% spanwise averaged buoyancy field


clear;%clc; 
addpath(genpath('/glade/u/home/liuchihl/TOOLS'));
% input coefficients
for cc=[1]
if cc==1
%1) D=0.5, Ri=0.16
  fname1 = 'doubleshearlayer/Ri016/R1/3D/D_0.5_';
  RI=0.16; Re=1000; Pr=1;
  LX=36.96; NX=768;
  LY=30;    NY=613;
  LZ=9.24;  NZ=192;
num=1;
%timestep = 400*10;
%  timestep = 400*17;
timestep = 400*15;
%timestep = 400*[2:24];

elseif cc==2
%2) D=1, Ri=0.16
fname1 = 'doubleshearlayer/Ri016/R1/3D/D_1_';
RI=0.16; Re=1000; Pr=1;
LX=78.54; NX=1536;
LY=30;    NY=613;
LZ=19.64;  NZ=384;
num=3;
timestep = 400*47;
%timestep = 400*[45:47];
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

%2) load the time of the 3d slices
base_dir=['/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/'];
%timestep=[4000:400:10000];
for i=1:num
   tt{i}=zeros(size(timestep));
   B_z{i}=zeros(NX,NY);
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
    filename_mean=[base_dir , fname1, num2str(i) '/mean.h5'];
    info = h5info([filename,timename]);
    %time
    tt{i}(k) = info.Groups.Attributes.Value;
    varname=['/Timestep/TH1'];

  if exist([filename,timename])
   B = h5read([filename,timename],varname);
   B_z{i} = mean(B,3);
  end
	end
end
varname=['/gyf/' '0001'];             % Y-COORDINATE
gyf(:)=h5read(filename_mean,varname);

% NX=size(U,1); NZ=size(U,3); NY=size(U,2);
x=linspace(0,LX,NX); y=gyf;%linspace(-LY(i)/2,LY(i)/2,NY); 
z=linspace(0,LZ,NZ);

if cc==1
%    save B_yavg_0.5.mat B_z y x tt
elseif cc==2
    save B_yavg_1.mat B_z y x tt 
elseif cc==3
    save B_yavg_2.mat B_z y x tt 
else cc==4
    save B_yavg_3.mat B_z y x tt 
end

end
