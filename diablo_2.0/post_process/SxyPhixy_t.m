% This script is part of the Butterfly effect project: it calculates the cross-spectra between upper and lower disturbances without smoothing the speectrum, here, unlike SxyPhixy.m, our goal is to calculate the phase difference with respect to time via 2D slices from the upper and lower layer
% see crossspec.m for squared coherence (smoothed)
% In this section, data is being extracted and variables are calculated and
% saved
clear;%clc;

% Define file names and directory
%fname1 = 'Ri_0.16_0.5_20_small/';
 base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_10/butterfly/Ri_0.14_0.05_';
% Input coefficients in DIABLO coordinate
RI=0.14; Re=1000;Pr=1;
LX=28.08; NX=512;
LY=20;    NY=361;
LZ=7.02;  NZ=128;
addpath(genpath('/glade/u/home/liuchihl/TOOLS'));
num=15;
nk=501;

% grids
x=linspace(0,LX,NX);z=linspace(0,LZ,NZ);
dx=mean(diff(x));
nx=length(x); nz=length(z);


%% cross-spectra between upper and lower layer
for i = 1:num
filename_mean=[base_dir num2str(i) '/mean.h5'];
filename=[base_dir num2str(i) '/movie.h5'];

%file_info=h5info(filename_mean);
%att_info=file_info.Groups.Attributes;
%nk=att_info.Value
nk=501;
x=linspace(0,LX,NX);z=linspace(0,LZ,NZ);
dx=mean(diff(x));
nx=length(x); nz=length(z);
time{i}=zeros(1,nk);
Sxyk{i} = zeros(1,nk);
phixyk{i} = zeros(1,nk);
%N=nx; 
kN = 1/(2*dx); % Nyquist wavenumber

for k=1:nk-200

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
 time{i}(k)=h5read(filename_mean,varname);

%  varname=['/v_xz/' timename];
 varnamePlus = ['/vPlus_xz/' timename];
 varnameMinus = ['/vMinus_xz/' timename];

%  v_xz=h5read(filename,varname);
 vPlus = h5read(filename,varnamePlus);
 vMinus = h5read(filename,varnameMinus);

% v' (vertical velocity perturbation)
vPlus_prime = vPlus-mean(vPlus,[1,2]);
vMinus_prime = vMinus-mean(vMinus,[1,2]);

cc=2;   % a choice between taking the avg after fft or taking avg before fft
	% 1: standard way: taking avg after fft
	% 2: taking the average of the series then take fft

	if cc == 1
% compute the cross correlation with the written function: fft_cpsd
for j=1:nz
	[Sxx,Syy,Sxy,kj] = fft_cpsd(vMinus_prime(:,j),vPlus_prime(:,j),dx,'rec');

% compute phase phi
	Cxy = real(Sxy);             % Co-spectrum
	Qxy = -imag(Sxy);            % Quadrature specetrum
	phixy = atan2(-Qxy,Cxy)*360/(2*pi); % phase in degrees

	if j==1
		temp1 = zeros(size(kj))';
		temp2 = zeros(size(kj))';
	end
	temp1 = temp1 + Sxy;
	temp2 = temp2 + phixy;
end
[~,po] = min(abs(kj*2*pi-.44));
% take the avg in the z (spanwise direction)
Sxyk{i}(k) = temp1(po)./nz;
phixyk{i}(k) = temp2(po)./nz;
%save phasediff_16 Sxyk phixyk time

	else

%2) compute the cross correlation by taking the average of the series first

     [Sxx,Syy,Sxy,kj] = fft_cpsd(mean(vMinus_prime,2),mean(vPlus_prime,2),dx,'rec');

% compute phase phi
        Cxy = real(Sxy);             % Co-spectrum
        Qxy = -imag(Sxy);            % Quadrature specetrum
        phixy = atan2(-Qxy,Cxy)*360/(2*pi); % phase in degrees

[~,po] = min(abs(kj*2*pi-.44));
	Sxyk{i}(k) = Sxy(po);
	phixyk{i}(k) = phixy(po);
%save phasediff_16_avg Sxyk phixyk time
	end
end 
i
end
save phasediff_12_avg Sxyk phixyk time




%% phase between subharmonic and KH mode at z=0 from 2D slice
nk=501;
x=linspace(0,LX,NX);z=linspace(0,LZ,NZ);
dx=mean(diff(x));
nx=length(x); nz=length(z);

   
for i = 1:num
filename_mean=[base_dir num2str(i) '/mean.h5'];
filename=[base_dir num2str(i) '/movie.h5'];
time{i}=zeros(1,nk);
dphi{i} = zeros(1,nk);
%N=nx; 
kN = 1/(2*dx); % Nyquist wavenumber

for k=1:nk-100

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
 time{i}(k)=h5read(filename_mean,varname);

 varname=['/v_xz/' timename];
 v_xz=h5read(filename,varname);
 
% v' (vertical velocity perturbation)
vprime = v_xz-mean(v_xz,[1,2]);

% taking the spanwise average of the series then take fft
  vpm = mean(vprime,2);                   % series
        % determine the wavenumber
        if mod(nx,2)==0;                       % if N is even
            kj = 2*pi*(-nx/2:nx/2-1)/(nx*dx);         % frequencies
        else mod(nx,2)==1;                     % if N is odd
            kj = 2*pi(-(nx-1)/2:(nx-1)/2)/(nx*dx);   % frequencies
        end
        % fft
        V=zeros(size(kj));
        V = fftshift(fft(vpm,length(kj)))/nx;   % Fourier coeeficient
        % consider only the right half of the wavenumber
        kj = kj(nx/2+1:nx); 
        V = V(nx/2+1:nx);
% compute phase phi for both 2nd(subharmonic) and 3rd(KH) mode
   %1) take the ratio of V_KH and V_sub
   V32 = V(3)./V(2).^2;
        C = real(V32);               % real part
        Q = imag(V32);               % imaginery part
        dphi{i}(k) = atan2(Q,C);            % phase in radians

% [~,po] = min(abs(kj*2*pi-.44));
% 	Sxyk{i}(k) = Sxy(po);
% 	phixyk{i}(k) = dphi(po);
%save phasediff_16_avg Sxyk phixyk time

end
i
end
save dphi_sub_KH_16_211 dphi time



%% careful diagnostic of phase between subharmonic and KH mode at z=0 from 2D slice
nk=300;
x=linspace(0,LX,NX);z=linspace(0,LZ,NZ);
dx=mean(diff(x));
nx=length(x); nz=length(z);

   
i =7;
filename_mean=[base_dir num2str(i) '/mean.h5'];
filename=[base_dir num2str(i) '/movie.h5'];
% time=zeros(1,nk);
dphi = zeros(1,nk);
%N=nx; 
kN = 1/(2*dx); % Nyquist wavenumber

% k=190 is when t=195.85, soon after pairing begins
 for k=1:nk

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

 varname=['/v_xz/' timename];
 v_xz=h5read(filename,varname);
 
% v' (vertical velocity perturbation)
vprime = v_xz-mean(v_xz,[1,2]);

% taking the spanwise average of the series then take fft
  vpm(:,k) = mean(vprime,2);                   % series
        % determine the wavenumber
        if mod(nx,2)==0;                       % if N is even
            kj = 2*pi*(-nx/2:nx/2-1)/(nx*dx);         % frequencies
        else mod(nx,2)==1;                     % if N is odd
            kj = 2*pi(-(nx-1)/2:(nx-1)/2)/(nx*dx);   % frequencies
        end
        % fft
        V = zeros(size(kj));
        V = fftshift(fft(vpm(:,k),length(kj)))/nx;   % Fourier coeeficient
        % consider only the right half of the wavenumber
        kj = kj(nx/2+2:nx); 
        V = V(nx/2+2:nx);
% compute phase phi for both 1st(subharmonic) and 2nd(KH) mode
   %1) take the ratio of V_KH and V_sub
        C1 = real(V(1));               % real part
        Q1 = imag(V(1));               % imaginery part
        phi_sub(k) = atan2(Q1,C1)/pi;            % phase in radians
        C2 = real(V(2));               % real part
        Q2 = imag(V(2));               % imaginery part
        phi_kh(k) = atan2(Q2,C2)/pi;            % phase in radians
        
   V21 = V(2)./V(1);
        C = real(V21);               % real part
        Q = imag(V21);               % imaginery part
        dphi(k) = atan2(Q,C)/pi;            % phase in radians
   V211 = V(2)./(V(1).^2);
        C211 = real(V211);               % real part
        Q211 = imag(V211);               % imaginery part
        dphi211(k) = atan2(Q211,C211)/pi;            % phase in radians
 end
 %% plot phi_kh, phi_sub and phi_kh-2*phi_sub

 % take off the jump between pi and -pi: dphi and phixyk

    dd = diff(dphi211);
    po = find(abs(dd)>1.8);
    po(po<3) = [];
    dphi211(po)= nan;
    
      dd = diff(phi_sub);
    po = find(abs(dd)>1.8);
    po(po<3) = [];
    phi_sub(po)= nan;
    
      dd = diff(phi_kh);
    po = find(abs(dd)>1.8);
    po(po<3) = [];
    phi_kh(po)= nan;

% color=[84, 68, 187;213, 21, 21;255, 128, 0]/255;
colors = lines(3);
close all;
    hh=figure('position',[50 100 700 320]);
    t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
    nexttile
    plot(time,phi_kh,'linewidth',2.2,'color',colors(1,:)); hold on;
    plot(time,phi_sub,'linewidth',2.2,'color',colors(2,:)); 
    plot(time,dphi211,'linewidth',2.2,'color',colors(3,:)); 
    plot([0,300],[0.5,0.5],'k--','linewidth',1.2);
    plot([0 nk],[0 0],'k');
set(gca,'fontsize',13,'xminortick','on','yminortick','on','tickdir','out',...
    'linewidth',1.2);
xlabel('t','fontsize',15,'fontname','Lucida Bright','fontangle','italic'); 
ylabel('\phi/\pi','fontsize',15,'fontname','Lucida Bright','fontangle',...
    'italic');
legend('\phi_{KH}','\phi_{sub}','\phi^{KH}_{sub}','fontsize',12,...
    'location','northoutside','orientation','horizontal','box','off','fontname','Lucida Bright',...
    'fontangle','italic')

xlim([0 300]);
 
    print -dpng -r200 phiKHsub_dphi_14_7.png
 
 
 %%
 tsub = 200;    %overturning time
 tindx = find(abs(time-tsub)==min(abs(time-tsub)));
     close all;
figure('position',[50,100,700,600]);
t = tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
nexttile
 plot(x,vpm(:,tindx),'linewidth',2.2,'color','k'); hold on;
 plot([0 LX],[0 0],'k');
 xlim([0 LX]);
 title(['Ri=0.14, case #' num2str(i,'%.0f') ', t=' num2str(time(tindx),'%.1f')],'fontsize',12','fontweight','normal');
xlabel('x','fontsize',18);
ylabel('<w''>_y','fontsize',18);
set(gca,'fontsize',12,'xminortick','on','yminortick','on');

nexttile 
color=[84, 68, 187;213, 21, 21;255, 128, 0]/255;

    plot(time,phi_kh,'linewidth',2.2,'color',color(1,:)); hold on;
    plot(time,phi_sub,'linewidth',2.2,'color',color(2,:)); 
    plot(time,dphi,'linewidth',2.2,'color',color(3,:)); 
    plot(time,dphi211,'linewidth',2.2,'color','k'); 
    plot([time(tindx),time(tindx)],[-1 1],'k--');
    plot([0 nk],[0 0],'k');
xlabel('t','fontsize',18); ylabel('\phi/\pi','fontsize',18);
legend('\phi_{kh}','\phi_{sub}','\Delta\phi_{1}','\Delta\phi_{2}','fontsize',12,'location','northoutside','orientation','horizontal')
set(gca,'fontsize',12,'xminortick','on','yminortick','on');
xlim([0 time(nk)]);


% print -dpng phiKHsub_dphi_14_6.png





%% time series of phase difference phi
% Ri=0.12
cases = [2,5,6,8,9,12];             % all interesting cases
%  Ri=0.14
% cases = [1,2,3,9,10,12];             % all interesting cases
% Ri=0.16
cases = [1,3,5,7,9,14];  

c16=load('phasediff_16_avg.mat');
c14=load('phasediff_14_avg.mat');
c12=load('phasediff_12_avg.mat');



% color1=cbrewer('seq','Blues',num+2);
% color1=color1(3:end,:); %this is just because the first couple of colours in this colourmap are really light
color2=lines(length(cases));
close all;
hh=figure('position',[50 100 800 400]);
% plots
for i=1:length(cases)
bb(i)=plot(c16.time{cases(i)}(1:200),c16.phixyk{cases(i)}(1:200)/180,'-',...
    'linewidth',2,'color',color2(i,:));hold on;
plot(c16.time{cases(i)}(1),c16.phixyk{cases(i)}(1)/180,'o','color',color2(i,:),'linewidth',2,'markersize',10);
end
hold on; plot([0 180],[0 0],'k');
plot([0 180],[0.35 0.35],'k');
set(gca,'xminortick','on','yminortick','on',...
    'ticklength',[0.015 0.015],'fontsize',12,'tickdir','out')
xlabel('$t$','fontsize',16,'interpreter','latex');
ylabel('$\phi_{lu}/\pi$','fontsize',16,'interpreter','latex');
xlim([0 180]);ylim([-.15 .4]);
% ylim([-1 1])
% legend(bb,'2','5','6','8','9','12','fontsize',12,'interpreter','latex','location','south');
% legend(bb,'1','2','3','9','10','12','fontsize',12,'interpreter','latex','location','southeast');
legend(bb,'1','3','5','7','9','14','fontsize',12,'interpreter','latex','location','southeast');

% print -djpeg -r200 PhivsTime12.jpg


%% initial phase difference (phi_i) vs locking-on time

k16 = load('butterflystats_16.mat','tke','Time');
k14 = load('butterflystats_14.mat','tke','Time');
k12 = load('butterflystats_12.mat','tke','Time');

num=15;
% first, get the initial phases
for i=1:num
    iniphase12(i) = c12.phixyk{i}(1)/180;
    iniphase14(i) = c14.phixyk{i}(1)/180;
    iniphase16(i) = c16.phixyk{i}(1)/180;
end
% second, find the lock-in time
%1) first approach: find the time of minimum TKE
for i=1:num
[~,po16(i)] = min(mean(k16.tke{i}));
[~,po14(i)] = min(mean(k14.tke{i}));
[~,po12(i)] = min(mean(k12.tke{i}));
t_lock16(i) = c16.time{i}(po16(i));  % time that starts to lock-on
t_lock14(i) = c14.time{i}(po14(i));  % time that starts to lock-on
t_lock12(i) = c12.time{i}(po12(i));  % time that starts to lock-on

end




% %2) second approch:
% % calculate the derivative of phi
% phi_optimal = [66,62,58];
% diff16 = zeros(num,219);
% diff14 = zeros(num,219);
% diff12 = zeros(num,219);
% 
% for i=1:num
% diff16 = diff(c16.phixyk{i}(1:220));
% diff14 = diff(c14.phixyk{i}(1:220));
% diff12 = diff(c12.phixyk{i}(1:220));
% 
% for j = 1:length(diff16)
%      % find the point and 10 later timestep that are all <1 degree
%      % Also within +-5 of the optimal phase
%      if diff16(j:j+10) < 1 &  c16.phixyk{i}(j)-phi_optimal(1) < 5
% lockpos16(i) = j;  
% break
%     end
% end
% for j = 1:length(diff14)
%      % find the point and 10 later timestep that are all <1 degree
%      % Also within +-5 of the optimal phase
%      if diff14(j:j+10) < 1 &  c14.phixyk{i}(j)-phi_optimal(2) < 5
% lockpos14(i) = j;  
% break
%     end
% end
% for j = 1:length(diff12)
%      % find the point and 10 later timestep that are all <1 degree
%      % Also within +-5 of the optimal phase
%      if diff12(j:j+10) < 1 &  c12.phixyk{i}(j)-phi_optimal(3) < 5
% lockpos12(i) = j;  
% break
%     end
% end
% 
% t_lock16(i) = c16.time{i}(lockpos16(i));  % time that starts to lock-on
% t_lock14(i) = c14.time{i}(lockpos14(i));  % time that starts to lock-on
% t_lock12(i) = c12.time{i}(lockpos12(i));  % time that starts to lock-on
% 
% end


%%
close all;
figure('position',[50,100,700,300]); 

axes('position',[.1,.2,.27,.7]);
    plot(iniphase16,t_lock16,'xk','linewidth',1.8); 
    set(gca,'xminortick','on','yminortick','on',...
    'ticklength',[0.015 0.015],'fontsize',12,'tickdir','out','xtick',[-1:.5:1]);
hold on; plot([-.15 -.15],[20 55],'k--');
plot([.35 .35],[20 55],'k--')
xlabel('$\Delta\phi_{i}/\pi$','fontsize',16,'interpreter','latex');
ylabel('$t_{lock}$','fontsize',16,'interpreter','latex');
% axis([-.22 .12 20 55]);
axis([-1 1 20 55])
title('$Ri_0=0.16$','interpreter','latex','fontsize',12);

axes('position',[.4,.2,.27,.7]);
    plot(iniphase14,t_lock14,'xk','linewidth',1.8); 
     set(gca,'xminortick','on','yminortick','on',...
    'ticklength',[0.015 0.015],'fontsize',12,'tickdir','out',...
    'yticklabel',[],'xtick',[-1:.5:1]);
hold on; plot([-.15 -.15],[20 55],'k--');
plot([.35 .35],[20 55],'k--')
xlabel('$\Delta\phi_{i}/\pi$','fontsize',16,'interpreter','latex');
% axis([-.22 .12 20 55]);
axis([-1 1 20 55])
title('$Ri_0=0.14$','interpreter','latex','fontsize',12);

axes('position',[.7,.2,.27,.7]);
    plot(iniphase12,t_lock12,'xk','linewidth',1.8); 
    set(gca,'xminortick','on','yminortick','on',...
    'ticklength',[0.015 0.015],'fontsize',12,'tickdir','out',...
    'yticklabel',[],'xtick',[-1:.5:1]);
hold on; plot([-.15 -.15],[20 55],'k--');
plot([.35 .35],[20 55],'k--');
xlabel('$\Delta\phi_{i}/\pi$','fontsize',16,'interpreter','latex');
%   axis([-.22 .12 20 55]);
axis([-1 1 20 55])
    title('$Ri_0=0.12$','interpreter','latex','fontsize',12);
print -djpeg -r200 phasediff_tlock_avg.jpg
% ylim([-5 5]);
% ylim([-.1 .1])
% xlim([20 90])










