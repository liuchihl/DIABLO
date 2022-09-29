% This script calculates the cross-spectra between upper and lower disturbances
% Here we can either use 2D slices or 3D flow field
%clear;%clc;

% input coefficients
fname1 = 'Ri_0.16_0.5_20_small/';
% In DIABLO coordinate
RI=0.16; Re=1000;Pr=1;
% d=10, LY=20, Ri=0.16, KICK=0.5
LX=28.28; NX=512;
LY=20;    NY=361;
LZ=7.03;  NZ=128;

timestep=[2000:2000:20000];
addpath(genpath('/glade/u/home/liuchihl/TOOLS'));
cm_1 = load('MPL_gnuplot.rgb');
i = 1;
U=zeros(NX,NY(i),NZ);V=zeros(NX,NY(i),NZ);
W=zeros(NX,NY(i),NZ);TH1=zeros(NX,NY(i),NZ);
v=zeros(NX,NY(i),NZ);

if     i==1
    base_dir=['/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_10/',fname1,'output/'];
end
x=linspace(0,LX(i),NX); y=linspace(-LY(i)/2,LY(i)/2,NY(i));
z=linspace(0,LZ(i),NZ);
dx=mean(diff(x));
nx=length(x); ny=length(y); nz=length(z);
%N=nx; 
kN = 1/(2*dx); % Nyquist wavenumber

filename=[base_dir];

for k=2%:length(timestep)

  if (timestep(k)<10)
    timename=['out0000' int2str(timestep(k)) '.h5'];
  elseif (timestep(k)<100)
    timename=['out000' int2str(timestep(k)) '.h5'];
  elseif (timestep(k)<1000)
    timename=['out00' int2str(timestep(k)) '.h5'];
  elseif (timestep(k)<10000)
    timename=['out0' int2str(timestep(k)) '.h5'];
  else
    timename=['out' int2str(timestep(k)) '.h5'];
  end
%varname1=['/Timestep/U'];
varname2=['/Timestep/V'];
%varname3=['/Timestep/W'];
%varname4=['/Timestep/TH1'];

if exist([filename,timename])
%U(:,:,:)=h5read([filename,timename],varname1);
V(:,:,:)=h5read([filename,timename],varname2);
%W(:,:,:)=h5read([filename,timename],varname3);
%TH1(:,:,:)=h5read([filename,timename],varname4);

info = h5info([filename,timename]);
%time
tt(k) = info.Groups.Attributes.Value;

% v' (vertical velocity perturbation)
v=V-mean(V,[1,3]);

% find the upper and lower layer that need to be crosscorrelated
[~,up]=min(abs(y-0.6585));
[~,low]=min(abs(y+0.6585));
vup = squeeze(v(:,up,:)); vlow = squeeze(v(:,low,:));
N=length(vup);
% compute the cross correlation with the written function: fft_cpsd
clear Coh phixy Cxy Qxy S_smoothxx S_smoothyy S_smoothxy kj_sm dphi phi_low phi_high sig
for j=1:nz
[Sxx,Syy,Sxy,kj] = fft_cpsd(vlow(:,j),vup(:,j),dx,'rec');

% here we smooth the PSD by band average with 6 dof
DOF=8;                         % Degree of freedom
M=DOF/2;                            % the length that is band averaged over with
for i=1:floor(length(kj)/M)
    S_smoothxx(i) = mean(Sxx(M*i-(M-1):M*i));
    S_smoothyy(i) = mean(Syy(M*i-(M-1):M*i));
    S_smoothxy(i) = mean(Sxy(M*i-(M-1):M*i));
         kj_sm(i) = mean([M*i-(M-1):M*i]/(N*dx));  
end
Coh = 4*conj(S_smoothxy).*S_smoothxy ./ (S_smoothxx.*S_smoothyy); % Smoothed squared Coherence

% confidence interval for coherence
Pro = .95; alpha = 1-Pro;   
qF=finv(Pro,2,2*M-2);                % this is in MATLAB definition
Co_cr= qF / ( (2*M-2)/2+ qF );       % critical value of coherence

% compute phase phi
Cxy = real(S_smoothxy);             % Co-spectrum
Qxy = -imag(S_smoothxy);            % Quadrature specetrum
phixy = atan2(-Qxy,Cxy)*360/(2*pi); % phase in degrees

% sig = find(Coh>=Co_cr)             % 95% significance
% confidence interval for the phase 
dphi = asin(sqrt( (2/(2*M-2).*((1-Coh)./Coh).*qF) ))*360/(2*pi);  % delta phi in degrees
phi_low = phixy-dphi; phi_high = phixy+dphi; 
if j==1
temp1 = zeros(1,floor(length(kj)/M));temp2 = zeros(1,floor(length(kj)/M));
temp3 = zeros(1,floor(length(kj)/M));temp4 = zeros(1,floor(length(kj)/M));
temp5 = zeros(size(kj))';
end
temp1 = temp1 + Coh;
temp2 = temp2 + phixy;
temp3 = temp3 + phi_low;
temp4 = temp4 + phi_high; 
temp5 = temp5 + Sxy;
end
% take the avg in the z (spanwise direction)
Coh = temp1./nz;
phixy = temp2./nz;
phi_low = temp3./nz;
phi_high = temp4./nz;
sig = find(Coh>=Co_cr);             % 95% significance
Sxy = temp5./nz;

% plots
%close all
hh=figure('position',[50 50 700 800]);
axes('position',[0.12 0.73 0.4 0.22]);
plot(x,mean(vup,2),'-k','linewidth',1.2); xlim([0 x(end)]);
hold on;plot(x,mean(vlow,2),'-r','linewidth',1.2); xlim([0 x(end)]);
xlabel('$X$','fontsize',12,'fontname','times','interpreter','latex'); 
ylabel({'$\left<w''\right>_y$'},'fontsize',12,'fontname','times','interpreter','latex'); 
set(gca,'fontsize',12,'xminortick','on','yminortick','on',...
    'ticklength',[0.015 0.015],'fontsize',12);
hold on; plot([0 x(end)],[0 0],'-k');
legend('Upper','Lower','fontsize',9,'fontname','times','location','southwest');
title(['$t$=' num2str(tt(k))],'fontweight','normal',...
    'interpreter','latex');


axes('position',[0.57 0.73 0.4 0.22]);
semilogy(tii,tke_int,'-k','linewidth',1.2); 
hold on; semilogy(tii(tii==tt(k)),tke_int(tii==tt(k)),'or' )
title('$TKE$','fontname','times','interpreter','latex');
xlabel('$t$','interpreter','latex','fontsize',12);
set(gca,'fontsize',12,'xminortick','on','yminortick','on',...
    'ticklength',[0.015 0.015],'fontsize',12);
xlim([0 tii(end)]);

axes('position',[0.12 0.41 0.8 0.25]);
aa=plot(kj_sm.*(2*pi),Coh,'k-','linewidth',1.1);
xlabel('$k^\ast$','interpreter','latex','fontsize',12);
ylabel('Coherence Squared $\bar{\hat{\gamma}}_{lu}^{2}(k^\ast)$','interpreter','latex','fontsize',12);
hold on; 
plot([0,kN.*(2*pi)],[Co_cr,Co_cr],'--k');
hold on; plot([kj(1)*2*pi,kj(1).*2*pi],[0,1],'k--');
%hold on; plot([1/6,1/6],[0,1],'k--');
text(kN.*2*pi*1.02,Co_cr,'95%','fontsize',12);
text((kN-.2)*(2*pi),.8,{'DOF=8'},'HorizontalAlignment', 'right',...
    'interpreter','latex','fontsize',12); 
axis([0 kN*(2*pi) 0 1]);
set(gca, 'YScale','linear','XScale','linear','xminortick','on','yminortick','on',...
    'ticklength',[0.015 0.015],'fontsize',12)

axes3=axes('position',[0.12 0.09 0.8 0.25]);
% plot(fj_sm(sig),phixy(sig),'k.','markersize',11);
for i=sig
errorbar(kj_sm(i)*(2*pi),phixy(i),abs(phixy(i)-phi_low(i)),...
    abs(phixy(i)-phi_high(i)),0,0,'.k','linewidth',1.0,'markersize',12);
hold on;
% text(2e-4*1.2,10.^0*S_smooth_tr(pos),'95%','fontsize',12,'color','r');
end
%hold on; plot([1/12,1/12],[-180,180],'k--');
%hold on; plot([1/6,1/6],[-180,180],'k--');
plot([0 kN*(2*pi)],[0 0],'-k','linewidth',.8);
plot([0 kN*(2*pi)],[0.65*pi*180/pi 0.65*pi*180/pi],'-k','linewidth',.8);
axis([0 kN*(2*pi) -180 180]);
xlabel('$k^\ast$','interpreter','latex','fontsize',12);
ylabel('Phase $\bar{\hat{\phi}}_{lu}(k^\ast)$ (deg)','interpreter','latex','fontsize',12);
set(gca,'ytick',[-180:120:180], 'YScale','linear','XScale','linear');
set(gca,'xminortick','on','yminortick','on',...
    'ticklength',[0.015 0.015],'fontsize',12)
axes3.YAxis.MinorTickValues = -180:30:180;
%saveas(hh,sprintf('CrossSpec%d.jpg',k));

end
end


%%% plot the cross spectrum (Sxy) instead of coherence square, because we want to see the data without smoothing, especially in the low wavenumbers

figure; 
%ax1 = axes('position',[0.1,0.1,.8,.8]);
semilogy(kj*2*pi,real(Sxy));

Cxy = real(Sxy);             % Co-spectrum
Qxy = -imag(Sxy);            % Quadrature specetrum
phixy = atan2(-Qxy,Cxy)*360/(2*pi); % phase in degrees
figure('position',[50,50,1000,500]);
ax1 = axes('position',[0.1,0.1,.8,.8]);
plot(kj*2*pi,phixy,'k','linewidth',1.2);hold on;
plot([0 kN*(2*pi)],[0 0],'-k','linewidth',.8);
plot([0 kN*(2*pi)],[0.65*pi*180/pi 0.65*pi*180/pi],'-k','linewidth',.8);
xlabel('$k^\ast$','interpreter','latex','fontsize',12);
ylabel('Phase $\bar{\hat{\phi}}_{lu}(k^\ast)$ (deg)','interpreter','latex','fontsize',12);
set(gca,'ytick',[-180:120:180], 'YScale','linear','XScale','linear');
set(gca,'xminortick','on','yminortick','on',...
    'ticklength',[0.015 0.015],'fontsize',12)
ax1.YAxis.MinorTickValues = -180:30:180;
ylim([-180,180]);
xlim([0 10]);






