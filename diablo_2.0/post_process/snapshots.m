% 2D snapshots:


addpath(genpath('/glade/u/home/liuchihl/TOOLS'));
% load('/glade/u/home/liuchihl/work/colormaps/storm_color.mat');
cm1 = load('NCV_rainbow2.rgb');% cm2 = load('/glade/u/home/liuchihl/work/colormaps/NCL_colormap/MPL_s3pcpn_l.rgb');
cm1 = cm1(46:211,:);

tt = [139,235,327,501];
close all;
% tt=1:501;
for i=1:3
    if i==1
    base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_10/butterfly/Ri_0.16_0.05_6';
    k=291;
    NY=361; LY=20;
    NX=512; LX=27.76;
    NZ=128; LZ=6.94;
    x=linspace(0,LX,NX); y=gyf;
    elseif i==2
    base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_10/butterfly/Ri_0.12_0.05_6';       
    k=221;
    NY=361; LY=20;
    NX=512; LX=28.28;
    NZ=128; LZ=7.07;
    x=linspace(0,LX,NX); y=gyf;
    else
    base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_10/butterfly/Ri_0.12_0.05_11';       
    k=177;
    NY=361; LY=20;
    NX=512; LX=28.28;
    NZ=128; LZ=7.07;
    x=linspace(0,LX,NX); y=gyf;
    end
filename=[base_dir '/mean.h5'];
    varname=['/gyf/' '0001'];             % Y-COORDINATE
    Z=h5read(filename,varname);
    X=x;
filename=[base_dir '/movie.h5'];

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
 varname1=['/th1_xy/' timename];
 Bxz=h5read(filename,varname1);

 if i==1
     figure('Position',[100,100,900,200]);
    t = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
        nexttile;
        pcolor(x,y,Bxz');shading interp;
        ylim([-6,6]); %xlabel('x','fontsize',18,'fontname','times','fontangle','italic');
        set(gca,'yticklabel',[],'xticklabel',[])
%         cmocean('thermal')
        caxis([-.8,.8])
 elseif i==2
         nexttile;
        pcolor(x,y,Bxz');shading interp;
        ylim([-6,6]);%xlabel('x','fontsize',18,'fontname','times','fontangle','italic');
        set(gca,'yticklabel',[],'xticklabel',[])
caxis([-.8,.8]) 
 else
         nexttile;
        pcolor(x,y,Bxz');shading interp;
        ylim([-6,6]);%xlabel('x','fontsize',18,'fontname','times','fontangle','italic');
        set(gca,'yticklabel',[])
        col = colorbar; set(col,'ytick',[-1:1:1],'yticklabel',[]);
        set(gca,'yticklabel',[],'xticklabel',[])
caxis([-.8,.8]) 
 end
end
 colormap(cm1/255)
% colormap jet
 print -dpng -r300 3snapshots.png
 
 
 %% 2D 2 snapshots
 
 addpath(genpath('/glade/u/home/liuchihl/TOOLS'));
% load('/glade/u/home/liuchihl/work/colormaps/storm_color.mat');
cm1 = load('NCV_rainbow2.rgb');% cm2 = load('/glade/u/home/liuchihl/work/colormaps/NCL_colormap/MPL_s3pcpn_l.rgb');
cm1 = cm1(46:211,:);
% tt = [232,227];
close all;
% tt=1:501;
for i=1:2
    if i==1
    base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_10/butterfly/Ri_0.14_0.05_12';
    k=231;
    NY=361; LY=20;
    NX=512; LX=28.08;
    NZ=128; LZ=7.02;
    x=linspace(0,LX,NX); %y=gyf;
    else i==2
    base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_10/butterfly/Ri_0.14_0.05_2';       
    k=227;
    NY=361; LY=20;
    NX=512; LX=28.08;
    NZ=128; LZ=7.02;
    x=linspace(0,LX,NX); %y=gyf;
    
    end
filename=[base_dir '/mean.h5'];
    varname=['/gyf/' '0001'];             % Y-COORDINATE
    Z=h5read(filename,varname);
    X=x; y=gyf;
filename=[base_dir '/movie.h5'];

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
 varname1=['/th1_xy/' timename];
 Bxz=h5read(filename,varname1);

 if i==1
     figure('Position',[100,100,700,200]);
axes('position',[.05,.05,.42,.9])
        pcolor(x,y,Bxz');shading interp;
        ylim([-4,4]); %xlabel('x','fontsize',18,'fontname','times','fontangle','italic');
        set(gca,'yticklabel',[],'xticklabel',[])
%         cmocean('thermal')
% caxis([-1.25,1.25])
caxis([-.8,.8])

 elseif i==2
        axes('position',[.51,.05,.42,.9])

        pcolor(x,y,Bxz');shading interp;
        ylim([-4,4]);%xlabel('x','fontsize',18,'fontname','times','fontangle','italic');
        set(gca,'yticklabel',[],'xticklabel',[])
%         cmocean('thermal')
co = colorbar('position',[.945,.05,.03,.9]); set(co,'ytick',[-1,1])
caxis([-.8,.8])
 end
end

 colormap(cm1/255)
% colormap jet
 print -dpng -r300 2snapshots.png
 