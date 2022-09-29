% run this code after readmean_h5.m
clear all;%close all;
addpath(genpath('/glade/u/home/liuchihl/TOOLS'));

% Ri=0.12
% fname1 = 'd_10/Ri_0.12';fname2 = 'd_5/Ri_0.12';
% fname3 = 'd_3/Ri_0.12/noslip';
% d = [14,5,3]; % the distance from the shear layer to the boundary
% NX=[1008,1024,1024];
% LX = [28.2,28.8,31.4];
% for j=1:3%:3
%     x{j}=linspace(0,LX(j),NX(j));
% end
% % delx = x{j}(2)-x{j}(1);
% nk=[1585,1565,1539];
% t_pick{1}=[211:20:411,416:20:604,623:20:705,707:20:905,907:20:1107,1109:1200];
% % t_pick{1}=[950];
% t_pick{2}=[200:20:360,362:20:611,613:20:862,864:20:1113];
% t_pick{3}=[200:20:251,253:20:502,504:20:753,755:20:1004,1006:20:1100];

% t_pick{1}=[1:200,211:411,413:604,623:705,707:905,907:1107,1109:1283,1285:1434,1436:nk(1)];
% t_pick{2}=[1:159,169:360,362:611,613:862,864:1113,1115:1364,1366:nk(2)];
% t_pick{3}=[1:251,253:502,504:753,755:1004,1006:1255,1257:nk(3)];

% Ri=0.16
fname1 = 'd_10/Ri_0.16';fname2 = 'd_5/Ri_0.16';
fname3 = 'd_3/Ri_0.16/';
d = [10,5,3]; % the distance from the shear layer to the boundary
NX=[1024,1024,1024];
LX = [28.2,27.9,28.8];
for j=1:3%:3
    x{j}=linspace(0,LX(j),NX(j));
end
nk=[1530,1403,1503];       %the final timestep
t_pick{1}=[100:20:301,303:20:627,629:20:953,955:20:1200];
t_pick{2}=[100:20:451,453:20:927,929:20:1200];
t_pick{3}=[100:20:301,303:20:902,904:20:1200];


% t_pick{1}=[1:301,303:627,629:953,955:1279,1281:nk(1)];
% t_pick{2}=[1:451,453:927,929:nk(2)];
% t_pick{3}=[1:301,303:902,904:nk(3)];

clear lambda_mean
for  j = 1:3%:3
    if j==1
base_dir=['/glade/campaign/univ/uosc0024/diablo_2.0/Large_case/',fname1];
    elseif j==2    
base_dir=['/glade/campaign/univ/uosc0024/diablo_2.0/Large_case/',fname2];
    else
base_dir=['/glade/campaign/univ/uosc0024/diablo_2.0/Large_case/',fname3];
    end
filename=[base_dir '/movie.h5'];
filename_mean=[base_dir '/mean.h5'];
% lsa_L = [14.2,14.4,15.7];
% init=300; fin = 800;
% clear lambda_mean h aspect WW E
for k=1:length(t_pick{j})
    if (t_pick{j}(k)<10)
        timename=['000' int2str(t_pick{j}(k))];
    elseif (t_pick{j}(k)<100)
        timename=['00' int2str(t_pick{j}(k))];
    elseif (t_pick{j}(k)<1000)
        timename=['0' int2str(t_pick{j}(k))];
    else
        timename=[int2str(t_pick{j}(k))];
    end

    varname=['/time/' timename];            % TIME
    time{j}(k)=h5read(filename_mean,varname); 

    varname=['/gyf/' timename];             % Y-COORDINATE
    gyf{j}=h5read(filename_mean,varname); 
    
    varname=['/th1_xy/' timename];
    
    TH1=h5read(filename,varname);
    TH1_y= mmderiv(gyf{j}',TH1')';%buoyancy (scalar) gradient
    
    varname=['/u_xy/' timename];
    U1 = h5read(filename,varname);
%do spectral analysis from -10 to +10, and average them out
% d1 = find(abs(gyf+10)==min(abs(gyf+10))); %the position at -10
% d2 = find(abs(gyf-10)==min(abs(gyf-10))); %the position at +10
% clear WW E lambda 
% A1_m=mean(TH1,2);                          % average over depth
% A1_y_m=mean(TH1_y,2);
Sp_s=0;clear Sp
% middle = round(length(gyf{j})/2);
% [~,upper] = min(abs(gyf{j}+gyf{j}(middle)-3)); 
% [~,lower] = min(abs(gyf{j}+gyf{j}(middle)+3));
% zz = lower:upper;

zz = 1:length(gyf{j});
for i=1:length(zz)
        yn = [TH1(:,zz(i));TH1(:,zz(i));TH1(:,zz(i));TH1(:,zz(i))];
        % the data
[Sp,kj]=fft_psd(yn,mean(diff(x{j})),'rec');
    Sp_s = Sp_s+Sp;                          % add all PSD together
end
Sp = Sp_s./i;
% [~,popo] = sort(Sp,'descend');

% figure; semilogx(kj,Sp);
%     ind = find(kj>1./x{j}(end)&kj<.5);%the wavenumber should be less than 0.5
    max_pos = find(Sp==max(Sp));
    lambda{j}(k)= 1./kj(max_pos);%wavelength
% [localmax,pomax] = findpeaks(A1_y_m); 
% [~,bb]= sort(localmax,'descend');
% maxposition = pomax(bb);
% if length(maxposition)>=2
%     L_est=abs(x{j}(maxposition)-x{j}(maxposition(1))); % distance from all the peaks to the largest peak
%     [~,po]=min(abs(L_est-lsa_L(j)));
%     lambda{j}(k) = L_est(po);
% 
% else   % if there's only 1 peak left
% lambda{j}(k) = abs(x{j}(maxposition(1))-x{j}(A1_y_m==min(A1_y_m)));
% end


% vertical scale of the billows
% first half
for xx=1:size(TH1,1)/2
    
    [~,inu] = min(abs((TH1(xx,:))-1.2));
    zU1{j}(xx,k) = gyf{j}(inu);
    [~,ind] = min(abs((TH1(xx,:))-0.8266));
    zL1{j}(xx,k) = gyf{j}(ind);
end
% vertical scale of the billows
% second half
xx=size(TH1,1)/2+1:size(TH1,1);
for i=1:length(xx)
    [~,inu] = min(abs((TH1(xx(i),:))-1.2));
    zU2{j}(xx(i)-xx(1)+1,k) = gyf{j}(inu);
    [~,ind] = min(abs((TH1(xx(i),:))-0.8266));
    zL2{j}(xx(i)-xx(1)+1,k) = gyf{j}(ind);
end

 
end
h{j} = max([max(zU1{j})-min(zL1{j});max(zU2{j})-min(zL2{j})]);

aspect{j}= h{j}./lambda{j};
end
%% plots
color=[84, 68, 187;213, 21, 21;255, 128, 0]/255;

figure;
set(gcf,'units','centimeters','paperunits','centimeters')
set(gcf,'PaperType','A4');
pp=[0.63 0.9 16 18];
ps=[0 0 pp(3)/1.1 pp(4)/1.1];
set(gcf,'paperposition',pp);
set(gcf,'position',ps);
cases = 1:3;

for j=cases
if j==1
    ax1 = axes('position',[.1 .76 .8 .2]);
end

plot(time{j},lambda{j},'-','linewidth',1.2,'color',color(j,:));
hold on;
set(gca,'fontsize',12,'xminortick','on','yminortick',...
    'on','xticklabel',[],'ticklength',[.015 .015]);
ylabel('$\lambda$','fontsize',12,'interpreter','latex');
% ylim([10 30]);xlim([time{3}(1) 280]);
ylim([10 30]);xlim([50 245]);
box on;
end
for j=cases
    if j==1
ax2 = axes('position',[.1 .54 .8 .2]);
    end
plot(time{j},h{j},'linewidth',1.2,'color',color(j,:));hold on;
set(gca,'fontsize',12,'xminortick','on','yminortick',...
    'on','xticklabel',[],'ticklength',[.015 .015]);
box on;
% ylim([0 5]);xlim([time{3}(1) 280]);
ylim([0 5]);xlim([50 245]);

ylabel('$h$','fontsize',12,'interpreter','latex');

end
for j = cases
if j==1
    ax3 = axes('position',[.1 .32 .8 .2]);
end
plot(time{j},aspect{j},'linewidth',1.2,'color',color(j,:));hold on;
set(gca,'fontsize',12,'xminortick','on','yminortick',...
    'on','xticklabel',[],'ticklength',[.015 .015]);
box on;
ylabel('$h/\lambda$','fontsize',12,'interpreter','latex');
% xlim([time{3}(1) 280]);
xlim([50 245]);

end
for j=cases
    if j==1
        ax4 = axes('position',[.1 .1 .8 .2]);
    end
plot(time{j},d(j)./lambda{j},'linewidth',1.2,'color',color(j,:));hold on;
set(gca,'fontsize',12,'xminortick','on','yminortick',...
    'on','ticklength',[.015 .015]);
box on;
ylabel('$d/\lambda$','fontsize',12,'interpreter','latex');
xlabel('time','fontsize',12);
% xlim([time{3}(1) 280]);
xlim([50 245]);
ylim([0 1]);
end
legend('$d$ = 10','$d$ = 5','$d$ = 3','interpreter','latex','location',...
    'northwest');
set(gcf,'color','w');
print -djpeg aspect_0.16.jpg


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
