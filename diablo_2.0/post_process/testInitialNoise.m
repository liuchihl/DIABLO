N_TH=1; % The number of scalars
Re = 1000; NU=1/Re; % Enter the Reynolds number or viscosity from input.dat
Pr=1;   kappa=NU/Pr; % Prandtl number
RI(1:N_TH)=0.16; % Enter the richardson number for each scalar

% In DIABLO coordinate 

%- (Ri=0.16, LY=20)
%NY=721; LY=20;
%NX=1024; LX=28.61;
%NZ=256; LZ=7.15;

%NY(2)=721; LY(2)=20;
%NX(2)=1024; LX(2)=28.61;
%NZ(2)=256; LZ(2)=7.15;

%- (Ri=0.16, LY=20, coarser)
NY=[361,361]; LY=20;
NX=[512,512]; LX=27.76;
NZ=[128,128]; LZ=6.94;


%- (Ri=0.16, LY=13.16)
%NX(2)=[1024]; NZ(2)=[256];  NY(2)=[397];         % Here, NY should match the value in grid_def.all
%LX(2)=[28.8]; LZ(2)=[7.2];  LY(2) = [13.16];     %vertical length

% NX=[128,128]; NY=[65,65]; NZ=[16 16];
% LX=[30,30]; LY=[20,20]; LZ=[4 4];

ys=[-7,-7];
% after determining the grid sizes and domain size, compute x,y,z

%    base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_3/Ri_0.16_0.5_20';
base_dir1='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_10/butterfly/tests/Ri_0.16_0.05_3';

%    base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_3/Ri_0.16_0.5_20_dup';
base_dir2='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_10/butterfly/tests/Ri_0.16_0.05_1';
% 2d slices
filename1=[base_dir1 '/movie.h5'];
filename2=[base_dir2 '/movie.h5'];
filename1mean=[base_dir1 '/mean.h5'];
filename2mean=[base_dir2 '/mean.h5'];

x=linspace(0,LX,NX(1)); y=linspace(-LY/2,LY/2,NY(1));
z=linspace(0,LZ,NZ(1));


% plot_xy = VideoWriter('Ri_.16,d=3,LY=20,KICK=0.5,compare2.mp4');
% 
% plot_xy.FrameRate = 20; plot_xy.Quality = 100;
% open(plot_xy);
close all
figure('position',[50,50,800,600])
for k=1%:2:400
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
   
varname1=['/u_xy/' timename];
U1=h5read(filename1,varname1); 
varname2=['/u_xy/' timename];
U2=h5read(filename2,varname2); 

varname1=['/th1_xy/' timename];
TH1=h5read(filename1,varname1); 
varname2=['/th1_xy/' timename];
TH2=h5read(filename2,varname2); 
    varname=['/time/' timename];            % TIME
    time1(k)=h5read(filename1mean,varname); 

    time2(k)=h5read(filename2mean,varname); 
u_p{1}=U1-mean(U1);
u_p{2}=U2-mean(U2);

%[~,po(i)]=min(abs(y-ys(i)));

% plots
 axes('position',[.1,.55,.8,.4]);
 pcolor(x,y,TH1'); shading flat
 caxis([0 2]);
%  axis equal tight
%  colormap(ax1,jet);
%  cmocean('thermal');
 % tit = title(sprintf('Time=%.2f',time(k)));
 set(gca,'tickdir','out','xminortick','on','yminortick','on','fontsize',12,...
     'xticklabel',[]);
 z1 = colorbar('southoutside'); set(z1,'position',[.06,.06,0.3,.015]);
 text(2,5,'$B$','fontsize',18,'interpreter','latex');
 ylabel('$Z$','fontsize',12,'interpreter','latex');
 title(sprintf('$Time=%.2f$',time1(k)),'interpreter','latex')

 axes('position',[.1,.1,.8,.4]);
pcolor(x,y,TH2'); shading flat
 caxis([0 2]);
%  axis equal tight
 colormap(jet);
%  cmocean('thermal');
 % tit = title(sprintf('Time=%.2f',time(k)));
 set(gca,'tickdir','out','xminortick','on','yminortick','on','fontsize',12);
 z1 = colorbar('southoutside'); set(z1,'position',[.06,.06,0.3,.015]);
 text(2,5,'$B$','fontsize',18,'interpreter','latex');
 ylabel('$Z$','fontsize',12,'interpreter','latex');
 xlabel('$X$','fontsize',12,'interpreter','latex');
%   title(sprintf('$Time=%.2f$',time2(k)),'interpreter','latex')

% drawnow
% M=getframe(gcf);
%  writeVideo(plot_xy,M);
% 
%  
% clf;
end
% close(plot_xy)


%%

close all;
% x vs rms(u')
figure('position',[20,150,1300,400]); 
for i=1:2
x=linspace(0,LX(i),NX(i)); y=linspace(-LY(i)/2,LY(i)/2,NY(i)); 
z=linspace(0,LZ(i),NZ(i));

plot(x,rms(u_p{i},2)); hold on;
xlabel('x','fontsize',16);title('$\sqrt{\left<u''^2\right>_y}(t=0)$','interpreter','latex','fontsize',18);
legend('LY=20','LY=13.16');
end
% y vs rms(u')
figure('position',[20,150,1300,400]); 
for i=1:2
x=linspace(0,LX(i),NX(i)); y=linspace(-LY(i)/2,LY(i)/2,NY(i)); 

plot(y(2:end),rms(u_p{i}(:,2:end),1)); hold on;
xlabel('y','fontsize',16);title('$\sqrt{\left<u''^2\right>_x}(t=0)$','interpreter','latex','fontsize',18);
legend('LY=20','LY=13.16');
end

% histogram
%close all
figure;
for i=1:2

% x=linspace(0,LX(i),NX(i)); y=linspace(-LY(i)/2,LY(i)/2,NY(i)); 

h=histogram(u_p{i}(:),60,'Normalization','pdf');     % counts
%h.FaceAlpha=0.2;
%p=histcounts(u_p{i}(:),60,'Normalization','pdf');
title('$u''(t=0)$','interpreter','latex','fontsize',18);
ylabel('PDF','fontsize',18);
% legend('LY=20','LY=13.16');
bin = (h.BinEdges(2:end)+h.BinEdges(1:end-1))/2;    % bin values, middle point of the bar

%u_p_av(i) = trapz(bin,bin.*pdf);
hold on;
grid minor;
set(gca,'fontsize',12','tickdir','out','xminortick','on','yminortick','on');
end
legend('LY=20','LY=13.16');
h.FaceAlpha=0.5;
mean(u_p{1}(:))
std(u_p{1}(:))

mean(u_p{2}(:))
std(u_p{2}(:))

%%


%%
% x vs u'
figure('position',[20,150,1300,400]); 
for i=1:2
plot(x,u_p{i}(:,po(i))); hold on;
xlabel('x');title('$u'', (t=0$, at shear layer)','interpreter','latex');
legend('LY=20','LY=13.16');
end
%print -djpeg pdf_u_p.jpg
% plot histogram of the u' value.

% figure; plot(y,rms(U-mean(U),1))
% 
% figure; pcolor(x,y,(U-mean(U))'); shading flat; 



