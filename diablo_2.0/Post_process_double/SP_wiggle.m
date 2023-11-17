% this script breaks down how tke shear production influenced by flow geometry
% we first extract the  3D flow field and calculate u', w' and u'w' at which the upper and lower shear layers are in phase and out of phase.
clear;
% D=2, Ri=0.16
% fname1 = 'doubleshearlayer/Ri016/R1/3D/D_2_';
% RI=0.16; Re=1000; Pr=1;
% LX=29.92; NX=576;
% LY=30;    NY=613;
% LZ=7.48;  NZ=144;
% num=1; 

fname1 = 'doubleshearlayer/Ri016/R1/3D/D_1_';
RI=0.16; Re=1000; Pr=1;
% LX=78.54; NX=1536;
% LY=30;    NY=613;
% LZ=19.64;  NZ=384;

num=1;
% timestep = 400*[49];
base_dir=['/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/'];

for i=1
filename=[base_dir,fname1, num2str(i) '/movie.h5'];
filename_mean=[base_dir,fname1, num2str(i) '/mean.h5'];
file_info=h5info(filename_mean);
att_info=file_info.Groups.Attributes;
nk=att_info.Value;
% k1 = 200:400;
k1 = 1:350;
x=linspace(0,LX,NX);
U=zeros(NX,NY);V=zeros(NX,NY);W=zeros(NX,NY);
u=zeros(NX,NY);v=zeros(NX,NY);w=zeros(NX,NY);
ume = zeros(NY,length(k1)); vme = zeros(NY,length(k1));

        for k=1:length(k1)

  if (k1(k)<10)
    timename=['000' int2str(k1(k))];
  elseif (k1(k)<100)
    timename=['00' int2str(k1(k))];
  elseif (k1(k)<1000)
    timename=['0' int2str(k1(k))];
  else
    timename=[int2str(k(k))];
  end
varname=['/time/' timename];            % TIME
time(k)=h5read(filename_mean,varname);
varname=['/gyf/' timename];             % vertical grids
gyf=h5read(filename_mean,varname);
varname=['/ume/' timename];             
ume(:,k)=h5read(filename_mean,varname);
varname=['/vme/' timename];             
vme(:,k)=h5read(filename_mean,varname);
       end

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
u(:,:,k)=U-repmat(ume(:,k)',[NX,1]);

varname=['/v_xy/' timename];
V=h5read(filename,varname);
v(:,:,k)=V-repmat(vme(:,k)',[NX,1]);

varname=['/th1_xy/' timename];
TH(:,:,k)=h5read(filename,varname);
        end
end

% out of phase
%t_out = 47;
% in phase
%t_in = 56;

%u_p = u(:,:,[t_out,t_in]);
%v_p = v(:,:,[t_out,t_in]);
%b = TH(:,:,[t_out,t_in]);
%ume = ume(:,[t_out,t_in]);
x = linspace(0,LX,NX);
y = gyf;

%save SP_wiggles.mat x y time u_p v_p b ume
% save SP_wiggles_time_3.mat x y time u v TH ume
% save('SP_wiggles_time_1_1.mat','x','y','time','u','v','TH','ume','-v7.3');

%figure('position',[1000,1000,800,700]); 
%t = tiledlayout(2,3,"TileSpacing","compact","Padding","compact");
%nexttile;
%pcolor(x,y,u(:,:,t_in)); shading interp;

%nexttile;
%pcolor(x,y,v(:,:,t_in)); shading interp;

%nexttile;
%pcolor(x,y,u(:,:,t_in).*v(:,:,t_in)); shading interp;


%nexttile;
%pcolor(x,y,v(:,:,t_in); shading interp;
%nexttile;
%pcolor(x,y,u(:,:,t_in); shading interp;
%nexttile;
%pcolor(x,y,u(:,:,t_in); shading interp;


%% movie of the dB/dz
% close all;
% load SP_wiggles_time_1_1.mat;
%%
close all
cm = load('NCV_rainbow2.rgb');
cm = cm(20:236,:);
plot_xy = VideoWriter('B_z_1_1.avi');
plot_xy.FrameRate = 20; plot_xy.Quality = 100;
open(plot_xy);

 figure('position',[50 50 600 750]);
 for i=20:300
Bz = mmderiv(y,TH(:,:,i)');
t = tiledlayout(1,4,"TileSpacing","compact","Padding","compact");
pcolor(x,y,Bz); shading interp;
hold on;
contour(x,y,Bz,[0.5:.1:1],'color','w','linewidth',1.2);
ylim([-2,2])
clim([0.5,1]);
colormap(cm/255);
set(gca,'tickdir','out','xminortick','on','yminortick','on','fontsize',12);
colorbar('southoutside');
% z1 = colorbar('southoutside'); set(z1,'position',[.06,.06,0.3,.015]);
text(2,1,'$dB/dz$','fontsize',18,'interpreter','latex','color','y');
ylabel('$Z$','fontsize',12,'interpreter','latex');
%ylim([-10 10 ]);
% ylim([0,6]);
xlabel('X','fontsize',12);%ylabel('Z','fontsize',12);
D=1;
title(['D=' num2str(D) ', Time=' num2str(time(i),'%.1f')],'fontname','times','fontweight','normal')
drawnow
M=getframe(gcf);
writeVideo(plot_xy,M); 
clf;
%  close
 end

close(plot_xy)
