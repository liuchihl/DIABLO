% This script computes the dissipation of tke (epsilon) and tpe (Chi) in 3D 


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
load('stats_16_0.5D.mat','thme','ume','eta_avg','Time','Rig');
[~,ind] = max(eta_avg{1}(1:500)); ind = ind+4;
%timestep = 400*10;
  timestep = 400*17;
%timestep = 400*[2:24];

elseif cc==2
%2) D=1, Ri=0.16
fname1 = 'doubleshearlayer/Ri016/R1/3D/D_1_';
RI=0.16; Re=1000; Pr=1;
LX=78.54; NX=1536;
LY=30;    NY=613;
LZ=19.64;  NZ=384;
num=1;
load('stats_16_D1.mat','thme','ume','eta_avg','Time');
[~,ind] = max(eta_avg{1}(1:1000));  ind = ind-4;
timestep = 400*49;
%timestep = 400*[30:50];
elseif cc==3
%3) D=2, Ri=0.16
fname1 = 'doubleshearlayer/Ri016/R1/3D/D_2_';
RI=0.16; Re=1000; Pr=1;
LX=29.92; NX=576;
LY=30;    NY=613;
LZ=7.48;  NZ=144; 
num=1;
load('stats_16_D2.mat','thme','ume','eta_avg','Time');
[~,ind] = max(eta_avg{1}(1:500)); ind = ind+3;
timestep = 400*37;
elseif cc==4
%4) D=3, Ri=0.16
fname1 = 'doubleshearlayer/Ri016/R1/3D/D_3_';
RI=0.16; Re=1000; Pr=1;
LX=28.56; NX=576;
LY=30;    NY=613;
LZ=7.14;  NZ=144;
num=1;
load('stats_16_D3.mat','thme','ume','eta_avg','Time');
[~,ind] = max(eta_avg{1}(1:500));  ind = ind-2;
 timestep = 400*39;
% timestep = 400*[10:54];
end

%X = linspace(0,LX,NX); Z = linspace(-LY/2,LY/2,NY);
%timestep=[4800,5600];  % d=10, case=3, at max KH amplitude

%timestep=[14800];
cm_1 = load('MPL_gnuplot.rgb');
% i = 1;
% num=1;
for m=num
TH1=zeros(NX,NY,NZ);

tt=zeros(size(timestep));
base_dir=['/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/'];

filename=[base_dir, fname1, num2str(m) '/output/'];
filename_mean=[base_dir , fname1, num2str(m) '/mean.h5'];
filename_movie=[base_dir , fname1, num2str(m) '/movie.h5'];
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
B(:,:,:)=h5read([filename,timename],varname4);

info = h5info([filename,timename]);
%time
tt(k) = info.Groups.Attributes.Value;

varname=['/gyf/' '0001'];             % Y-COORDINATE
	gyf(:)=h5read(filename_mean,varname);

% NX=size(U,1); NZ=size(U,3); NY=size(U,2);
x=linspace(0,LX,NX); y=gyf;%linspace(-LY(i)/2,LY(i)/2,NY); 
z=linspace(0,LZ,NZ);
% mean profiles
Bme = thme{1}(:,ind);
%Bme = sort(Bme);		%1) sort the buoyancy profile
By = mmderiv(y,Bme);		
	%2) without sorting
By = mmderiv(y,Bme);
U = ume{1}(:,ind);
b_p = B-mean(B,[1,2]);          % perturbation buoyancy

% buoyancy gradients
 for j=1:NY
   b_x(:,j,:) = gradient(squeeze(b_p(:,j,:))',mean(diff(x)))';
 end
 for j=1:NX
  b_y(j,:,:) = gradient(squeeze(b_p(j,:,:))',mean(diff(y)))';
  b_z(j,:,:) = gradient(squeeze(b_p(j,:,:)),mean(diff(z)));
 end
end
b_var = b_x.^2+b_y.^2+b_z.^2;

Ri=0.16; Re=1000; Pr=1;
%Chi = Ri/(Re*Pr) * b_var ./ repmat(By',[NX,1,NZ]);
%Chi_z0 = Chi(:,:,1);
Chi = 2*Ri/(Re*Pr) * b_var;
Chi_z0 = Chi(:,:,1);
%By_rep = repmat(By',[768,1,192]);
%ind1 = find(abs(y-1.5)==min(abs(y-1.5)));
%ind2= find(abs(y+1.5)==min(abs(y+1.5)));
%Chi = Ri/(Re*Pr) * b_var ./ mean(By_rep(:,ind2:ind1,1),[1,2]);

% extract epsilon from movie.h5
k=ind;
 if (k<10)
     timename=['000' int2str(k)];
   elseif (k<100)
     timename=['00' int2str(k)];
   elseif (k<1000)
     timename=['0' int2str(k)];
   else
     timename=[int2str(k)];
   end
varname2=['/epsilon_xy/' timename];
EPS=h5read(filename_movie,varname2);

end
end

save mixing_chi_eps_0.5_1_Bvariance.mat EPS Chi_z0 x y tt 
