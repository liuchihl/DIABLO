% save x-t data for all d cases
addpath(genpath('/glade/u/home/liuchihl/TOOLS'));
d=[10,4,3,2.5,2];

for cc=5%1:length(d)
	if cc==1
%1) 0.12_10
RI=0.12; Re=1000;Pr=1;
LX=28.28; NX=512;
LY=20;    NY=361;
LZ=7.07;  NZ=128;
base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_10/boundary/Ri_0.12_0.1_2';

	elseif cc==2
%2) 0.12_4
 RI=0.12; Re=1000;Pr=1;
 LX=28.21; NX=512;
 LY=20;    NY=361;
 LZ=7.05;  NZ=128;
 base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_4/boundary/Ri_0.12_0.1_3';

%3) 0.12_3
	elseif cc==3
RI=0.12; Re=1000;Pr=1;
LX=28.36; NX=512;
LY=20;    NY=361;
LZ=7.09;  NZ=128;
base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_3/boundary/Ri_0.12_0.1_2';

%4) 0.12_2.5
	elseif cc==4
RI=0.12; Re=1000;Pr=1;
LX=27.93; NX=512;
LY=20;    NY=361;
LZ=6.98;  NZ=128;
base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_2.5/Ri_0.12_0.1_2';

%5) 0.12_2    
    else
RI=0.12; Re=1000;Pr=1;
LX=29.16; NX=512;
LY=20;    NY=361;
LZ=7.29;  NZ=128;
base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_2/Ri_0.12_0.1_2';
    end
filename = [base_dir '/movie.h5'];
filename_mean=[base_dir '/mean.h5'];
varname=['/gyf/' '0001'];             % Y-COORDINATE
gyf=h5read(filename_mean,varname);

x=linspace(0,LX,NX); y=gyf; z=linspace(0,LZ,NZ);

[~,cen(cc)] = min(abs(y-d(cc)+LY/2));

 for k=1:501
 k
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
    varname1=['/u_xy/' timename];
    U=h5read(filename,varname1);
    Uxt(:,k) = U(:,cen(cc));
 end
 if cc==1
    save HovU10 Uxt x y time
 elseif cc==2
    save HovU4 Uxt x y time
 elseif cc==3
    save HovU3 Uxt x y time
 elseif cc==4
    save HovU25 Uxt x y time
 else cc==5
    save HovU2 Uxt x y time
 end

end
