% Run after readmean.m
% this xz means the horizontal plane (x and y) 
LX=28.36;
NX=512;
LZ=7.09;
NZ=128;
addpath(genpath('/glade/u/home/liuchihl/TOOLS'));

filename=[base_dir '/movie.h5'];

x=linspace(0,LX,NX);
z=linspace(0,LZ,NZ);
% plot_xz = VideoWriter('Butterfly_Ri_12_2uv_xz.avi');
plot_xz = VideoWriter('Bw_xy3_1.avi');

plot_xz.FrameRate = 20; plot_xz.Quality = 100;
open(plot_xz);
% close all;
% figure;
% set(gcf,'units','centimeters','paperunits','centimeters')
% set(gcf,'PaperType','A4');
% pp=[0.63 0.9 19.7 19.7];
% ps=[0 0 pp(3)/1.1 pp(4)/1.1];
% set(gcf,'paperposition',pp);
% set(gcf,'position',ps);

%%
close all;
    figure('position',[50,50,700,700]);

for k=1:200
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

% varname=['/th1_xz/' timename];
% A_th=h5read(filename,varname);
% for j=1:size(A_th,2)
%    A_th(:,j)=A_th(:,j)+drhodz*z(j);
%  end
% calculate u'w' (Reynolds stress in x-dir and z-dir)
varname=['/u_xz/' timename];
Uxz=h5read(filename,varname);
Umean=mean(Uxz,[1,2]); u=Uxz-Umean;

varname=['/v_xz/' timename];
Vxz=h5read(filename,varname);
Vmean=mean(Vxz,[1,2]); v=Vxz-Vmean;

uv=u.*v;
axes('position',[.1,.48,.8,.48])
pcolor(x,z,uv'); shading interp;
xlabel('$x$','interpreter','latex'); ylabel('$y$','interpreter','latex'); 
set(gca,'FontName','Times','FontSize',14,'ticklength',[.02,.02],...
    'xminortick','on','yminortick','on','tickdir','out');

title(['$u''w''$, $t$=' num2str(tii(k))],'fontweight','normal',...
    'interpreter','latex');

% caxis([-1e-3 0]);
colormap(jet(256));
% cmocean('deep')
colorbar;
hold on; contour(x,z,uv',[0 0],'k');
axes('position',[.2,.085,.6,.3])
 semilogy(tii,tke_int,'k'); hold on;
 semilogy(tii(k),tke_int(k),'ro');
xlabel('$t$','interpreter','latex'); ylabel('TKE','interpreter','latex'); 
set(gca,'FontName','Times','FontSize',14,'ticklength',[.02,.02],...
    'xminortick','on','yminortick','on','tickdir','out');

drawnow
M=getframe(gcf);
 writeVideo(plot_xz,M);
clf;

end

close(plot_xz);


%%
close all;
    figure('position',[50,50,900,500]);

for k=50:350
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

varname=['/th1_xz/' timename];
th=h5read(filename,varname);

varname=['/v_xz/' timename];
Vxz=h5read(filename,varname);

t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
nexttile;
pcolor(x,z,th'); shading interp;
xlabel('$x$','interpreter','latex'); ylabel('$y$','interpreter','latex'); 
set(gca,'FontName','Times','FontSize',16,'ticklength',[.02,.02],...
    'xminortick','on','yminortick','on','tickdir','out','linewidth',1.2);
title(['$B$, $t$=' num2str(tii(k))],'fontweight','normal',...
    'interpreter','latex');
box on
caxis([-1 1]);
cmocean('thermal')
colorbar;

t2 = nexttile;
pcolor(x,z,Vxz'); shading interp;
xlabel('$x$','interpreter','latex'); %ylabel('$y$','interpreter','latex'); 
set(gca,'FontName','Times','FontSize',16,'ticklength',[.02,.02],...
    'xminortick','on','yminortick','on','tickdir','out','linewidth',1.2);
title('$w$','fontweight','normal',...
    'interpreter','latex');
caxis([-.5 .5]);
colormap(t2,'jet')
colorbar

drawnow
M=getframe(gcf);
 writeVideo(plot_xz,M);
clf;

end

close(plot_xz);




%% movie of u'w' of the x-y (horizontal slice at z=0)
close all;
    figure('position',[50,50,800,600]);
   

tt=[5,13,23,49];

t = tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');

for k=1:length(tt)
  if (tt(k)<10)
    timename=['000' int2str(tt(k))];
  elseif (tt(k)<100) 
    timename=['00' int2str(tt(k))];
  elseif (tt(k)<1000)
    timename=['0' int2str(tt(k))];
  else
    timename=[int2str(tt(k))];
  end
varname=['/u_xz/' timename];
Uxz=h5read(filename,varname);
Umean=mean(Uxz,[1,2]); u=Uxz-Umean;

varname=['/v_xz/' timename];
Vxz=h5read(filename,varname);
Vmean=mean(Vxz,[1,2]); v=Vxz-Vmean;

uv=u.*v;
if k==1
    axes('position',[.1,.55,.36,.4]);
    pcolor(x,z,uv'); shading interp;
    set(gca,'FontName','Times','FontSize',14,'ticklength',[.02,.02],...
    'xminortick','on','yminortick','on','tickdir','out',...
    'xticklabel',[]);
    title(['$t$=' num2str(time(tt(k)),'%.2f')],'fontweight','normal',...
    'interpreter','latex');
    caxis(max(abs(caxis))*[-1 1])
    cb=colorbar('position',[.47,.55,.015,.4])
    set(cb,'Yaxis','right')
    text(.8,6.6,'(a)','fontsize',12)
elseif k==2
    axes('position',[.55,.55,.36,.4]);
    pcolor(x,z,uv'); shading interp;
    set(gca,'FontName','Times','FontSize',14,'ticklength',[.02,.02],...
    'xminortick','on','yminortick','on','tickdir','out',...
    'xticklabel',[],'yticklabel',[]);
    title(['$t$=' num2str(time(tt(k)),'%.2f')],'fontweight','normal',...
    'interpreter','latex');
    caxis(max(abs(caxis))*[-1 1])
    cb=colorbar('position',[.92,.55,.015,.4])
    set(cb,'Yaxis','right');
    text(.8,6.6,'(b)','fontsize',12)
elseif k==3
    axes('position',[.1,.08,.36,.4]);
    pcolor(x,z,uv'); shading interp;
    set(gca,'FontName','Times','FontSize',14,'ticklength',[.02,.02],...
    'xminortick','on','yminortick','on','tickdir','out');
    title(['$t$=' num2str(time(tt(k)),'%.2f')],'fontweight','normal',...
    'interpreter','latex');
    caxis(max(abs(caxis))*[-1 1])
    cb=colorbar('position',[.47,.08,.015,.4])
    set(cb,'Yaxis','right');
    text(.8,6.6,'(c)','fontsize',12)
else
    axes('position',[.55,.08,.36,.4]);
    pcolor(x,z,uv'); shading interp;
set(gca,'FontName','Times','FontSize',14,'ticklength',[.02,.02],...
    'xminortick','on','yminortick','on','tickdir','out',...
    'yticklabel',[]);
    title(['$t$=' num2str(time(tt(k)),'%.2f')],'fontweight','normal',...
    'interpreter','latex');
    caxis(max(abs(caxis))*[-1 1])
    cb=colorbar('position',[.92,.08,.015,.4])
    set(cb,'Yaxis','right');
    text(.8,6.6,'(d)','fontsize',12);
%     xlabel('x','FontName','Times','FontSize',14,'fontangle','italic');
end


end
colormap jet
% hold on; contour(x,z,uv',[0 0],'k');


% print -djpeg -r300 uwplot_14_12.jpg



%% save u and w at x-z and x-y slices
LX=28.08;
NX=512;
LZ=7.02;
NZ=128;
LY=20;
NY=361;
addpath(genpath('/glade/u/home/liuchihl/TOOLS'));
base_dir=['/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_10/butterfly/Ri_0.14_0.05_12'];

filename_mean=[base_dir '/mean.h5'];
filename=[base_dir '/movie.h5'];

x=linspace(0,LX,NX);
y=linspace(-10,10,NY);
z=linspace(0,LZ,NZ);
nk=501;

tt=1:nk;  % all time
Uxy=zeros(NX,NZ,nk);
Wxy=zeros(NX,NZ,nk);
Uxz=zeros(NX,NY,nk);
Wxz=zeros(NX,NY,nk);


for k=1:length(tt)
  if (tt(k)<10)
    timename=['000' int2str(tt(k))];
  elseif (tt(k)<100) 
    timename=['00' int2str(tt(k))];
  elseif (tt(k)<1000)
    timename=['0' int2str(tt(k))];
  else
    timename=[int2str(tt(k))];
  end
  
     varname=['/time/' timename];            % TIME
     t(k)=h5read(filename_mean,varname); 

varname=['/u_xz/' timename];
Uxy(:,:,k)=h5read(filename,varname);

varname=['/v_xz/' timename];
Wxy(:,:,k)=h5read(filename,varname);

varname=['/u_xy/' timename];
Uxz(:,:,k)=h5read(filename,varname);

varname=['/v_xy/' timename];
Wxz(:,:,k)=h5read(filename,varname);


end
X=x; Y=z; Z=y;
save uw_14_12 Uxy Wxy Uxz Wxz X Y Z t