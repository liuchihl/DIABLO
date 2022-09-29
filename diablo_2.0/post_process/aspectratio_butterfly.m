% run this code after readmean_h5.m
clear all;%close all;
addpath(genpath('/glade/u/home/liuchihl/TOOLS'));

 NY=361; LY=20;
 NX=512; LX=28.28;
 Re = 1000; NU=1/Re; % Enter the Reynolds number or viscosity from input.dat
Pr=1;   kappa=NU/Pr; % Prandtl number
RI=0.12; % Enter the richardson number for each scalar
num=15; 		 % number of cases
nk=501;
for  j = 1:num
base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_10/butterfly/Ri_0.12_0.05_';
filename_mean=[base_dir num2str(j) '/mean.h5'];
filename_movie=[base_dir num2str(j) '/movie.h5'];

% lsa_L = [14.2,14.4,15.7];
% init=300; fin = 800;
% clear lambda_mean h aspect WW E
for k=1:400
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
    time{j}(k)=h5read(filename_mean,varname); 

    varname=['/gyf/' timename];             % Y-COORDINATE
    gyf=h5read(filename_mean,varname); 
    
    varname=['/th1_xy/' timename];
    TH1=h5read(filename_movie,varname);
    
    %TH1_y= mmderiv(gyf',TH1')';%buoyancy (scalar) gradient
    
    varname=['/u_xy/' timename];
    U1 = h5read(filename_movie,varname);
    x=linspace(0,LX,NX);
%do spectral analysis from -10 to +10, and average them out
% d1 = find(abs(gyf+10)==min(abs(gyf+10))); %the position at -10
% d2 = find(abs(gyf-10)==min(abs(gyf-10))); %the position at +10
% clear WW E lambda 
% A1_m=mean(TH1,2);                          % average over depth
% A1_y_m=mean(TH1_y,2);
Sp_s=0;clear Sp
middle = round(length(gyf)/2);
[~,upper] = min(abs(gyf+gyf(middle)-6)); 
[~,lower] = min(abs(gyf+gyf(middle)+6));
zz = lower:upper;

zz = 1:length(gyf);
for i=1:length(zz)
        yn = TH1(:,zz(i));
        % the data
[Sp,kj]=fft_psd(yn,mean(diff(x)),'rec');
    Sp_s = Sp_s+Sp;                          % add all PSD together
end
Sp = Sp_s./i;
% [~,popo] = sort(Sp,'descend');

% figure; semilogx(kj,Sp);
%     ind = find(kj>1./x{j}(end)&kj<.5);%the wavenumber should be less than 0.5
    [Sp_sort,po] = sort(Sp,'descend');
    kj_max = kj(po(1));
    lambda{j}(k)= 1./(kj_max);%wavelength


% vertical scale of the billows
% for xx=1:size(TH1,1)
%     [~,inu] = min(abs((TH1(xx,:))-0.29));
%     zU{j}(xx,k) = gyf(inu);
%     [~,ind] = min(abs((TH1(xx,:))+.29));
%     zL{j}(xx,k) = gyf(ind);
% end


% first half
for xx=1:size(TH1,1)/2
    
    [~,inu] = min(abs((TH1(xx,:))-0.29));
    zU1{j}(xx,k) = gyf(inu);
    [~,ind] = min(abs((TH1(xx,:))+.26));
    zL1{j}(xx,k) = gyf(ind);
end
% vertical scale of the billows
% second half
xx=size(TH1,1)/2+1:size(TH1,1);
for i=1:length(xx)
    [~,inu] = min(abs((TH1(xx(i),:))-0.29));
    zU2{j}(xx(i)-xx(1)+1,k) = gyf(inu);
    [~,ind] = min(abs((TH1(xx(i),:))+0.26));
    zL2{j}(xx(i)-xx(1)+1,k) = gyf(ind);
end

 
end
% h{j} = max(abs(zU{j}-zL{j}));
h{j} = max([max(zU1{j})-min(zL1{j});max(zU2{j})-min(zL2{j})]);

aspect{j}= h{j}./lambda{j};
end

save aspect12 aspect h lambda time
%% plots
color=[84, 68, 187;213, 21, 21;255, 128, 0]/255;

figure;
set(gcf,'units','centimeters','paperunits','centimeters')
set(gcf,'PaperType','A4');
pp=[0.63 0.9 16 18];
ps=[0 0 pp(3)/1.1 pp(4)/1.1];
set(gcf,'paperposition',pp);
set(gcf,'position',ps);
cases = 9;

for j=1:length(cases)
if j==1
    ax1 = axes('position',[.1 .76 .8 .2]);
end

plot(time{cases(j)},lambda{cases(j)},'-','linewidth',1.2,'color',color((j),:));
hold on;
set(gca,'fontsize',12,'xminortick','on','yminortick',...
    'on','xticklabel',[],'ticklength',[.015 .015]);
ylabel('$\lambda$','fontsize',12,'interpreter','latex');
% ylim([10 30]);xlim([time{3}(1) 280]);
ylim([10 30]);xlim([0 300]);
box on;
end
for j=length(cases)
    if j==1
ax2 = axes('position',[.1 .54 .8 .2]);
    end
plot(time{cases(j)},h{cases(j)},'linewidth',1.2,'color',color((j),:));hold on;
set(gca,'fontsize',12,'xminortick','on','yminortick',...
    'on','xticklabel',[],'ticklength',[.015 .015]);
box on;
% ylim([0 5]);xlim([time{3}(1) 280]);
% ylim([0 6]);
xlim([0 300]);

ylabel('$h$','fontsize',12,'interpreter','latex');

end
for j = length(cases)
if j==1
    ax3 = axes('position',[.1 .32 .8 .2]);
end
plot(time{cases(j)},aspect{cases(j)},'linewidth',1.2,'color',color((j),:));hold on;
set(gca,'fontsize',12,'xminortick','on','yminortick',...
    'on','ticklength',[.015 .015]);
box on;
ylabel('$h/\lambda$','fontsize',12,'interpreter','latex');
% xlim([time{3}(1) 280]);
xlim([0 300]);
ylim([0 .4])
end

% legend('$d$ = 10','$d$ = 5','$d$ = 3','interpreter','latex','location',...
%     'northwest');
set(gcf,'color','w');
% print -djpeg aspect_0.16.jpg


%% calculate the spanwise wavelength with spectral analysis
% LX=41.88;
% NX=384;
LY=20;
NY = 501;
LZ=4;
NZ=128;

filename=[base_dir '/movie.h5'];
% x=linspace(0,LX,NX);
z=linspace(0,LZ,NZ);
clear lambda_mean
lam = [];
for j = 1:3
    if j==1
base_dir=['/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/',fname1];
    elseif j==2    
base_dir=['/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/',fname2];
    else
base_dir=['/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/',fname3];
    end
filename=[base_dir '/movie.h5'];
pos = find(abs(gyf+d(j))<0.0001);
t1 = 160; t2 = 660;%t1:initial time to measure KH; t2: the end time
tt = t1:t2;
clear lambda_mean h aspect WW E
for k=1:length(tt)%217:nk
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
varname=['/th1_zy/' timename];

% varname=['/th1_xz/' timename];
clear lambda
A2=h5read(filename,varname);
%do spectral analysis for all y, and average them out
for l = 1:length(gyf)
    [WW,E]=ezfft(z,A2(:,l),'space'); %WW: wavenumber, E:power spectrum
%     ind = find(WW<.5&WW>.2);%the wavenumber should be less than 0.7
    max_pos = find(E.*WW.^2==max(E.*WW.^2));
    lambda(l) = 2*pi./WW(max_pos);%wavelength
end
lambda_mean(k) = mean(lambda);
end

lam = [lam;lambda_mean];
end
clear lambda_mean 
lambda_mean = lam;
figure; 
for j=1:3
% if j==1
%     ax1 = axes('position',[.1 .76 .8 .2]);
% end
plot(time(t1:t2),lambda_mean(j,:),'o-','markersize',1.5,'linewidth',2);hold on;
set(gca,'fontsize',12,'xminortick','on','yminortick','on');
ylabel('\lambda_y','fontsize',16);xlabel('Time','fontsize',12);
%ylim([10 22]);
xlim([time(t1) time(t2)]);
end
legend('d_b = 10','d_b = 5','d_b = 3');

set(gcf,'color','w');

%%
figure;

plot(time(t1:t2),h1(1,:),'b');
hold on;plot(time(t1:t2),-h2(1,:),'b-');
figure;
plot(time(t1:t2),h1(2,:)+5,'r');
hold on;plot(time(t1:t2),-h2(2,:)-5,'r-');
figure;
plot(time(t1:t2),h1(3,:)+7,'k');
hold on;plot(time(t1:t2),-h2(3,:)-7,'k-');
