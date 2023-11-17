 clear
addpath(genpath('/glade/u/home/liuchihl/TOOLS'));
load('/glade/u/home/liuchihl/work/colormaps/storm_color.mat');
% cm = load('NCV_jaisnd.rgb')./255;
cm = load('NCV_rainbow2.rgb');

% base_dir=['/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_2/Ri_0.12_0.1_1'];
% base_dir=['/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_10/butterfly/Ri_0.12_0.05_8'];

base_dir=['/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_3/boundary/Ri_0.12_0.1_3'];
filename_mean=[base_dir '/mean.h5'];
    varname=['/gyf/' '0001'];             % Y-COORDINATE
    gyf=h5read(filename_mean,varname);

%     Ri=0.12; d=10
% NY=361; LY=20;
% NX=512; LX=28.28;
% NZ=128; LZ=7.07;

%     Ri=0.12, d=3
    NY=361; LY=20;
    NX=512; LX=28.36;
    NZ=128; LZ=7.09;
i=1;
x=linspace(0,LX(i),NX(i)); y=gyf;%y=linspace(-LY(i)/2,LY(i)/2,NY(i)); 
z=linspace(0,LZ(i),NZ(i));
filename=[base_dir '/movie.h5'];
clear plot_xy M
plot_xy = VideoWriter('Boundary_3_3');

% plot_xy = VideoWriter('Boundary_0.12_2_2.avi');
% 
plot_xy.FrameRate = 28; plot_xy.Quality = 100;
open(plot_xy);


% moviename = 'B_moving_2_1.gif';
tt=20:210;
close all;
    h=figure('position',[50,50,600,500]);
t=tiledlayout(1,1,'TileSpacing','compact','Padding','compact');
    nexttile;
  for k=tt
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

varname1=['/th1_xy/' timename];

 
TH=h5read(filename,varname1);

pcolor(x,gyf,TH');shading interp
caxis([-1 1]);
%  axis equal tight
shading flat
z2 = colorbar; 
colormap(cm(30:256-29,:)./255);
%cmocean('ice')
text(2,9,'$B$','fontsize',18,'interpreter','latex','color','y');
xlabel('X','fontsize',18,'fontname','times','fontangle','italic');
ylabel('Z','fontsize',12,'fontname','times','fontangle','italic');
set(gca,'tickdir','out','xminortick','on','yminortick','on','fontsize',18,'linewidth',1.2);
title(sprintf(['d=3, t=%.2f'],time(k)),'fontname','times','fontangle','italic','fontweight','normal','fontsize',15)
set(gcf,'color','w');
 drawnow

% % Capture the plot as an image 
%       frame = getframe(h); 
%       im = frame2im(frame); 
%       [imind,cmm] = rgb2ind(im,256); 
%       % Write to the GIF File 
%       if k == tt(1) 
%           imwrite(imind,cmm,moviename,'gif', 'Loopcount',inf,'DelayTime',.2); 
%       else 
%           imwrite(imind,cmm,moviename,'gif','WriteMode','append','DelayTime',.08); 
%       end 


 drawnow
M=getframe(gcf);
 writeVideo(plot_xy,M); 
% clf;
%  close
 end
% 
 close(plot_xy)

%

 % create an animated gif instead of a avi video
