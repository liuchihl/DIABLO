% save buoyancy slices
clear
tt = [1,114,167,197,245,501];       % d=10
 tt=[1,139,203,268,501];	     % d=2.5
% tt=1:501
clear Bxz Byz
% tt=1:501;
%base_dir=['/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_10/boundary/Ri_0.12_0.1_3'];
base_dir=['/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_2.5/Ri_0.12_0.1_1'];
filename_mean=[base_dir '/mean.h5'];
    varname=['/gyf/' '0001'];             % Y-COORDINATE
%     gyf=h5read(filename_mean,varname);
%  LX=28.28; NX=512;
%  LY=20;    NY=361; 	 
%  LZ=7.07;  NZ=128;
  
   LX=27.93; NX=512;
   LY=20;    NY=361;
   LZ=6.98;  NZ=128;
    i=1;
x=linspace(0,LX(i),NX(i)); %y=gyf;%y=linspace(-LY(i)/2,LY(i)/2,NY(i)); 
z=linspace(0,LZ(i),NZ(i));
    Z=h5read(filename_mean,varname);
    X=x;
    Y=z;
filename=[base_dir '/movie.h5'];

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
 Bxz(:,:,k)=h5read(filename,varname1);
%  varname2=['/th1_zy/' timename];
%  Byz(:,:,k)=h5read(filename,varname2);
%  varname3=['/th1_xz/' timename];
%  Bxy(:,:,k)=h5read(filename,varname3);
 
    varname=['/time/' timename];            % TIME
    time(k)=h5read(filename_mean,varname)
 end
 t=time;
save B_xz_12_25_1 Bxz X Y Z t
