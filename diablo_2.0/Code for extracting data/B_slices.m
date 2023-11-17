% save buoyancy slice
base_dir=['/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/'];
clear Bxz Byz

for  cc=[2]
if cc==1
%1) D=0.5, Ri=0.16
  fname1 = 'doubleshearlayer/Ri016/R1/3D/D_0.5_';
  RI=0.16; Re=1000; Pr=1;
  LX=36.96; NX=768;
  LY=30;    NY=613;
  LZ=9.24;  NZ=192;
num=2;
%tt = [1,141,276,382,676];
tt = [146,192,293,536]; % for CASE #2
%tt=136;   % tKH
elseif cc==2
%2) D=1, Ri=0.16
fname1 = 'doubleshearlayer/Ri016/R1/3D/D_1_';
RI=0.16; Re=1000; Pr=1;
LX=78.54; NX=1536;
LY=30;    NY=613;
LZ=19.64;  NZ=384;
num=1;
tt = [443,492,628,938];
%tt=415;   %tKH for DLM
%tt = 224; %tKH for SLM
elseif cc==3
%3) D=1.5, Ri=0.16
fname1 = 'doubleshearlayer/Ri016/R1/3D/D_1.5_';
RI=0.16; Re=1000; Pr=1;
LX=31.416; NX=576;
LY=30;    NY=613;
LZ=7.854;  NZ=144;
num=1;
%tt = [1,362,440,504,611];
tt = 311;
elseif cc==4
%3) D=2, Ri=0.16
fname1 = 'doubleshearlayer/Ri016/R1/3D/D_2_';
RI=0.16; Re=1000; Pr=1;
LX=29.92; NX=576;
LY=30;    NY=613;
LZ=7.48;  NZ=144;
num=1;
tt = [362,440,504,611]
%tt=370;
elseif cc==5
%4) D=3, Ri=0.16
fname1 = 'doubleshearlayer/Ri016/R1/3D/D_3_';
RI=0.16; Re=1000; Pr=1;
LX=28.56; NX=576;
LY=30;    NY=613;
LZ=7.14;  NZ=144;
num=1;
%tt = [1,379,470,500,631];
tt=387;
end

filename_movie=[base_dir, fname1, num2str(num) '/movie.h5'];
filename_mean=[base_dir , fname1, num2str(num) '/mean.h5'];
varname=['/gyf/' '0001'];             % Y-COORDINATE
gyf(:)=h5read(filename_mean,varname);

x=linspace(0,LX,NX); y=gyf;%linspace(-LY(i)/2,LY(i)/2,NY); 
z=linspace(0,LZ,NZ);
Bxz = zeros(length(x),length(y),length(tt));
Byz = zeros(length(z),length(y),length(tt));
 for k=1:length(tt)
%  figure('position',[50 50 1300 750])

 k
   if (tt(k)<10)
     timename=['000' int2str(tt(k))];
   elseif (tt(k)<100) 
     timename=['00' int2str(tt(k))];
   elseif (tt(k)<1000)
     timename=['0' int2str(tt(k))];
   else
     timename=[int2str(tt(k))];
   end
   
 varname1=['/th1_xy/' timename];
 Bxz(:,:,k)=h5read(filename_movie,varname1);
 varname2=['/th1_zy/' timename];
 Byz(:,:,k)=h5read(filename_movie,varname2);
    varname=['/time/' timename];            % TIME
    time(k)=h5read(filename_mean,varname)
 end   
 t=time;
if cc==1
save B_xz_16_0.5_2.mat Bxz Byz x y z  t
elseif cc==2
save B_xz_16_1_1.mat Bxz Byz x y z t
elseif cc==3
save B_xz_16_1.5_1_tkh.mat Bxz Byz x y z t
elseif cc==4
save B_xz_16_2_1.mat Bxz Byz x y z t
else cc==5
save B_xz_16_3_1_tkh.mat Bxz Byz x y z t
end

end
