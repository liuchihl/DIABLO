%% this script is to calculate the strain rate and the rotation rate for which to analyze elliptical instability
clear ; 
addpath(genpath('/glade/u/home/liuchihl/TOOLS'));
load('/glade/u/home/liuchihl/work/colormaps/storm_color.mat');

base_dir=['/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_10/boundary/Ri_0.12_0.1_3'];
filename_mean=[base_dir '/mean.h5'];
    varname=['/gyf/' '0001'];             % Y-COORDINATE
    gyf=h5read(filename_mean,varname);

% Ri=0.12


% cases with LY=20 
% NY=361; LY=20;
% NX=512; LX=28.28;
% NZ=128; LZ=7.07;

% 
%  NY=361; LY=20;
%  NX=512; LX=28.21;
%  NZ=128; LZ=7.05;
% 
% NY=361; LY=20;
% NX=512; LX=28.36;
% NZ=128; LZ=7.09;
% 
% LX=27.93; NX=512;
% LY=20;    NY=361;
% LZ=6.98;  NZ=128;

LX=29.16; NX=512;
LY=20;    NY=361;
LZ=7.29;  NZ=128;

% test (small)
%fname1= 'KH_test2';
%LX=30; LY=20; LZ=4;
%NX=128;NY=65; NZ=16;
% timestep=[100:100:4000];
% after determining the grid sizes and domain size, compute x,y,z

x=linspace(0,LX,NX); y=gyf;
z=linspace(0,LZ,NZ);

%% save data
filename=[base_dir '/movie.h5'];
tt=1:400
for k=tt
%  figure('position',[50 50 1300 750])

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
    varnameU=['/u_xy/' timename];
    U=h5read(filename,varnameU);
    varnameV=['/v_xy/' timename];
    V=h5read(filename,varnameV);
    
    r12 = mmderiv(x,V)-mmderiv(y,U')';
    e11 = mmderiv(x,U); e22 = mmderiv(y,V')';
    e12 = 1/2*(mmderiv(y,U')'+mmderiv(x,V));
% compute eij*eij and rij*rij
Str(:,:,k) = e11.^2+e22.^2+e12.^2+e12.^2;  % strain rate tensor
Rot(:,:,k) = 1/4*(r12.^2+r12.^2);


end
save strain_rotation12_10_3 Str Rot time x y




%% plots 
figure('position',[100,100,800,500]);
col = lines(5);
for i=1:5 
    if i==1
    load  strain_rotation12_10_3
    elseif i==2
        load  strain_rotation12_4_3
    elseif i==3
        load  strain_rotation12_3_2
    elseif i==4
        load  strain_rotation12_25_3
    else 
        load  strain_rotation12_2_1
    end
    
    a(i)=plot(time,squeeze(mean(Str,[1,2])),'linewidth',2,'color',col(i,:)); hold on
    b(i)=plot(time,squeeze(mean(Rot,[1,2])),'--','linewidth',2,'color',col(i,:)); hold on
    
end
set(gca,'fontsize',16,'xminortick','on');
xlim([50 300]);
legend(a,'d=10','d=4','d=3','d=2.5','d=2','fontsize',12);
xlabel('t'); ylabel('e_{ij}e_{ij}, r_{ij}r_{ij}/4','fontangle','italic','fontname','times');
set(gcf,'color','w')
%% animation
%close all;
clear plot_xy M
filename=[base_dir '/movie.h5'];
% plot_xy = VideoWriter('Butterfly9');

plot_xy = VideoWriter('Boundary_0.12_2_2.avi');
% 
plot_xy.FrameRate = 25; plot_xy.Quality = 100;
open(plot_xy);

%%
tt= 50:460%400;%100:2:380;%:1100%:500%:280;%:500;
% tt= 1052:1382;
% for k=1:501
%     if (k<10)
%         timename=['000' int2str(k)];
%     elseif (k<100)
%         timename=['00' int2str(k)];
%     elseif (k<1000)
%         timename=['0' int2str(k)];
%     else
%         timename=[int2str(k)];
%     end
% 
%     varname=['/time/' timename];            % TIME
%     time(k)=h5read(filename_mean,varname); 
% end
close all;
    figure;
set(gcf,'units','centimeters','paperunits','centimeters')
set(gcf,'PaperType','A4');
pp=[0.63 1 38 19.8];
ps=[0 0 pp(3)/1.1 pp(4)/1.1];
set(gcf,'paperposition',pp);
set(gcf,'position',ps);
 for k=tt
%  figure('position',[50 50 1300 750])

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
% varname = ['/ume_z/' timename];
% UME_Z = h5read(filename,varname); 
varname1=['/th1_xy/' timename];
varname2=['/epsilon_xy/' timename];
varnameU=['/u_xy/' timename];
%  varname4=['/u_xz/' timename];
% varnamePlus = ['/uPlus_xz/' timename];
%  varnameMinus = ['/uMinus_xz/' timename];
 
TH=h5read(filename,varname1);
EPS=h5read(filename,varname2);
U=h5read(filename,varnameU);
% U_xz=h5read(filename,varname4); 
% uPlus = h5read(filename,varnamePlus); 
% uMinus = h5read(filename,varnameMinus); 
% A1=mmderiv(gyf,U')';
% A2=mmderiv(gyf,TH')';
A3 = EPS; 
% [A3(k),po(k)] = max(mean(EPS));
%  
%   
ax1 = axes('position',[.06 .15 .3 .7]);

pcolor(x,gyf,U'); shading flat
caxis([-1 1]);
%  axis equal tight
colormap(ax1,jet);
% tit = title(sprintf('Time=%.2f',time(k)));
set(gca,'tickdir','out','xminortick','on','yminortick','on','fontsize',12);
z1 = colorbar('southoutside'); set(z1,'position',[.06,.06,0.3,.015]);
text(2,9,'$U$','fontsize',18,'interpreter','latex');
ylabel('$Z$','fontsize',12,'interpreter','latex');
% ylim([-5 5 ]);

ax2 = axes('position',[.37 .15 .3 .7]);
pcolor(x,gyf,TH');shading interp
%  caxis([-1 4]);
% caxis([0 2]);
caxis([-1 1]);
%  axis equal tight
shading flat
set(gca,'yticklabel',[],'tickdir','out','xminortick','on','yminortick','on','fontsize',12);
z2 = colorbar('southoutside'); set(z2,'position',[.37,.06,0.3,.015]);
colormap(ax2,flipud(cm))
%cmocean('ice')
text(2,9,'$B$','fontsize',18,'interpreter','latex');
xlabel('X','fontsize',12);%ylabel('Z','fontsize',12);
title(sprintf(['$Pr=1, Ri_0=0.12, d=2$, case 2, $t=%.2f$'],time(k)),'interpreter','latex')

ax3 = axes('position',[.68 .15 .3 .7]);
pcolor(x,gyf,log10(A3)'); shading flat
caxis([-6 -2]);
%  axis equal tight;
colormap(ax3,cm_eps/255);
text(2,9,'$\epsilon$','fontsize',25,'color','y','interpreter','latex');
% ylim([-5 5 ]);
set(gca,'yticklabel',[],'tickdir','out','xminortick','on','yminortick','on','fontsize',12);
z3 = colorbar('southoutside'); set(z3,'position',[.68,.06,0.3,.015]);

clear M
%  if size(M.cdata,1)~=1306||size(M.cdata,2)~=680
% %      
% %  set(gcf,'units','centimeters','paperunits','centimeters')
% %  set(gcf,'PaperType','A4');
% %  pp=[0.63 1 38 19.8];
% %  ps=[0 0 pp(3)/1.1 pp(4)/1.1];
% %  set(gcf,'paperposition',pp)
% %  set(gcf,'position',ps)
%  drawnow
%   end
 drawnow
M=getframe(gcf);
 writeVideo(plot_xy,M); 
clf;
%  close
 end
% 
 close(plot_xy)