clear all;
% this xz means the horizontal plane (x and y)
%LX=27.73;
%NX=512;
%LZ=6.93;
%NZ=128;
% Input coefficients in DIABLO coordinate
RI=0.16; Re=1000;Pr=1;
LX=27.76; NX=512;
LY=20;    NY=361;
LZ=6.94;  NZ=128;
addpath(genpath('/glade/u/home/liuchihl/TOOLS'));
base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_10/butterfly/Ri_0.16_0.05_';
num=15;
% base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/KH_test2'

for j=1:num
j
filename=[base_dir num2str(j) '/movie.h5'];
filename_mean=[base_dir num2str(j) '/mean.h5'];
file_info=h5info(filename_mean);
att_info=file_info.Groups.Attributes;
nk=att_info.Value;

x=linspace(0,LX,NX);
z=linspace(0,LZ,NZ);
U=zeros(NX,NZ);V=zeros(NX,NZ);W=zeros(NX,NZ);
u=zeros(NX,NZ);v=zeros(NX,NZ);w=zeros(NX,NZ);

time{j}=zeros(1,nk); 
Sp_ku = zeros(2,nk);Sp_kv = zeros(2,nk);
Sp_kw = zeros(2,nk); TKE{j} = zeros(2,nk);
%
for k=1:nk
k;
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

varname=['/u_xz/' timename];
U=h5read(filename,varname);
Umean=mean(U,[1,2]); u=U-Umean;

varname=['/v_xz/' timename];
V=h5read(filename,varname);
Vmean=mean(V,[1,2]); v=V-Vmean;

varname=['/w_xz/' timename];
W=h5read(filename,varname);
Wmean=mean(W,[1,2]); w=W-Wmean;

%TKE = .5*(u.^2+v.^2+w.^2);                    % turbulent kinetic energy at y=0

% PSD 
dx = mean(diff(x)); kN = 2*pi/(2*dx);         % Nyquist wavenumber
Sp_su=0;Sp_sv=0;Sp_sw=0;
for i=1:size(U,2)                           % loop over z direction 
     [Sp_u,kj]=fft_psd(u(:,i),dx,'rec');
     Sp_su = Sp_su+Sp_u;                          % add all PSD together
     [Sp_v,kj]=fft_psd(v(:,i),dx,'rec');
     Sp_sv = Sp_sv+Sp_v;
     [Sp_w,kj]=fft_psd(w(:,i),dx,'rec');
     Sp_sw = Sp_sw+Sp_w;
end
% if we want the time series of the entire spectrum
% Sp_ku(:,k) = Sp_su'./NZ;                          % take average
% Sp_kv(:,k) = Sp_sv'./NZ;
% Sp_kw(:,k) = Sp_sw'./NZ;

% if we only want the time series of the 2 lowest wavenumbers
Sp_ku(:,k) = Sp_su(1:2,:)./NZ;                          % take average
Sp_kv(:,k) = Sp_sv(1:2,:)./NZ;
Sp_kw(:,k) = Sp_sw(1:2,:)./NZ;

end
TKE{j} = .5*(Sp_ku+Sp_kv+Sp_kw);
end

save TKEspecta_lowestwavenumber16 TKE time
%% plot spectrogram
%close all;
% figure('position',[50,50,800,500]); 
% axes('position',[.1,.12,.775,.85]);
% [cc,hh]=contourf(time,log10(kj*2*pi),log10(TKE),[-40:1:-10,-10:.01:-1]);set(hh,'edgecolor','none');
% hold on; plot([0 time(end)],[log10(0.44),log10(0.44)],'k--');
%[cc,hh]=contourf(time,(kj*2*pi),log10(Sp_k),[-40:1:-10,-10:.01:-1]);set(hh,'edgecolor','none');
%hold on; plot([0 time(end)],[(0.44),(0.44)],'k--');

%xlabel('$t$','fontsize',18,'interpreter','latex');
%ylabel('$\log(k^\ast)$','fontsize',18,'interpreter','latex');
%col = colorbar('position',[.9 .2 .015 .6]);
%ylabel(col,'TKE spectra','fontsize',12)
%cbarrow
%colormap(jet);
%caxis([-10,-1]);
%set(gca,'fontsize',12,'xminortick','on','yminortick','on','tickdir','out');
%grid minor
%ylim([log10(2*pi/LX),1.6])
%print -djpeg spectrogram_butterfly_2.jpg


%% TKE spectra at fastest growing mode (k*=0.44) and subharmonic mode (k*~0.23)
% FGM and subharmonic mode are the 2nd and 1st lowest wavenumber
num=15;
% Ri=0.12
load TKEspecta_lowestwavenumber12
% cases = [2,5,6,8,9,12];             % all interesting cases
%  Ri=0.14
% load TKEspecta_lowestwavenumber14
% cases = [1,2,3,9,10,12];             % all interesting cases
% Ri=0.16
% load TKEspecta_lowestwavenumber16
% cases = [1,3,5,7,9,14];             % all interesting cases
% cases = [1:15];

color1=cbrewer('seq','Blues',num+2);
color1=color1(4:end,:); %this is just because the first couple of colours in this colourmap are really light
color2=cbrewer('seq','Reds',num+2);
color2=color2(4:end,:);
close all;
f=figure('position',[50,50,500,410]);
axes('position',[.165,.565,.8,.41]);
for j=cases%1:num
a(j)=plot(time{j},TKE{j}(1,:),'linewidth',2);%,'color',color1(j,:));
hold on;
end
set(gca,'fontsize',12,'xticklabel',[],'xminortick','on','yminortick','on','tickdir','out');
ylabel('$E(k=0.22)$','fontsize',12,'interpreter','latex');
xlim([50 400]);
legend('1','3','5','7','9','14','fontsize',12,'interpreter','latex','location','northwest');

axes('position',[.165,.125,.8,.41]);
for j=cases%1:num
b(j)=plot(time{j},TKE{j}(2,:),'linewidth',2);%,'color',color2(j,:));
hold on;
end
set(gca,'fontsize',12,'xminortick','on','yminortick','on','tickdir','out');
xlabel('$t$','fontsize',12,'interpreter','latex');
ylabel('$E(k=0.44)$','fontsize',12,'interpreter','latex');

xlim([50 400]);

% print -djpeg -r200 TKE_spec_subharmonic_KHmode_16.jpg



%% Growth rate of the fastest growing mode (k*=0.44) and subharmonic mode (k*~0.23)
% with TKE spectra: sigma = 1/2 dln(E)/dt
num=15;
% Ri=0.12
% load TKEspecta_lowestwavenumber12
% cases = [2,5,6,8,9,12];             % all interesting cases
%  Ri=0.14
% load TKEspecta_lowestwavenumber14
% cases = [1,2,3,9,10,12];             % all interesting cases
% Ri=0.16
load TKEspecta_lowestwavenumber16
cases = [1,3,5,7,9,14];             % all interesting cases
for j=1:num
sig1{j} = 1/2*mmderiv(time{j},log10(TKE{j}(2,:)));
sig_sub{j} = 1/2*mmderiv(time{j},log10(TKE{j}(1,:)));
end

color1=cbrewer('seq','Blues',num+2);
color1=color1(4:end,:); %this is just because the first couple of colours in this colourmap are really light
color2=cbrewer('seq','Reds',num+2);
color2=color2(4:end,:);
close all;
f=figure('position',[50,50,500,410]);
axes('position',[.165,.565,.8,.41]);
for j=cases%1:num
a(j)=plot(time{j},sig_sub{j},'linewidth',2);%,'color',color1(j,:));
hold on;
end
set(gca,'fontsize',12,'xticklabel',[],'xminortick','on','yminortick','on','tickdir','out');
ylabel('$\sigma_{sub}$','fontsize',12,'interpreter','latex');
legend('1','3','5','7','9','14','fontsize',12,'interpreter','latex','location','northeast');
xlim([0 250]);ylim([-0.04 .05])


axes('position',[.165,.125,.8,.41]);
for j=cases%1:num
b(j)=plot(time{j},sig1{j},'linewidth',2);%,'color',color2(j,:));
hold on;
end
set(gca,'fontsize',12,'xminortick','on','yminortick','on','tickdir','out');
xlabel('$t$','fontsize',12,'interpreter','latex');
ylabel('$\sigma_{1}$','fontsize',12,'interpreter','latex');

xlim([0 250]);ylim([-0.04 .05])

% print -djpeg -r200 GrowthRateSubharmonic_mode1_16.jpg
%% note:
% Ri=12:
%  most of the cases pairs, 
%  6,12 doesn't pair but primary energy is comparable to other cases 
%  8 is interesting, where the energy is much smaller for the primary
%  instability, but has large subharmonic energy
%  9 has similar strength of primary energy but lower subharmonic energy,
%  result in weaker turbulence

% Ri=0.14:
% 3 weak pairing, 4 strong pairing
% 5, 7,8 pairs, 6 weak pairing
% 9, 12 pairs, 10 no pairing, 11 weak pairing
% 13 no pairing, 14 weak pairing, 15 strong pairing

%  Ri=0.16:
% 1 weak primary inst., 2-6 same energy of primary inst.
% 9,11 weak primary, 7,8,10,12 same energy of primary inst.,
% 13,15 same energy, 14 weaker primary inst. energy


%% PSD with -5/3
%colors = jet(length(1:10:501));
%close all;
%figure; 

%for i=1%:length(colors)
%aa(i)=loglog(2*pi*kj,TKE(:,10*i-9),'color',colors(i,:));hold on;
%end
%hold on; 
%plot(2*pi*(kj),0.000001*(kj).^(-5/3),'--k','linewidth',1.1);
%axis([1e-1 5e1 1e-10 1e0])
%set(gca,'fontsize',12);

%xlabel('k*');ylabel('PSD of TKE')

%figure;
%aa(i)=loglog(2*pi*kj,TKE(:,63),'o-','color','k');hold on;
%hold on; 
%plot(2*pi*(kj),0.000001*(kj).^(-5/3),'--k','linewidth',1.1);



