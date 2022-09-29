% This script calculates the cross-spectra between upper and lower disturbances without smoothing the speectrum
% Here we can either use 2D slices or 3D flow field
% see crossspec.m for squared coherence (smoothed)
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

for k=1:3

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

% the position of the upper and lower layer, which is symmetric about y(vertical)=0
h0 = linspace(0,3,100);
for hh=1:length(h0)
% find the upper and lower layer that need to be crosscorrelated
[~,up]=min(abs(y-h0(hh)));
[~,low]=min(abs(y+h0(hh)));
vup = squeeze(v(:,up,:)); vlow = squeeze(v(:,low,:));
N=length(vup);
% compute the cross correlation with the written function: fft_cpsd
for j=1:nz
	[Sxx,Syy,Sxy,kj] = fft_cpsd(vlow(:,j),vup(:,j),dx,'rec');

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
% take the avg in the z (spanwise direction)
%Sxyk{hh} = temp1./nz;
phixyk = temp2./nz;
[~,po] = min(abs(2*pi*kj-0.44)); % FGM
phixy_optimal{k}(:,hh) = phixyk(po);
end  %for height hh
end  %for time k
end  %for i cases

close all;
hh=figure('position',[50 50 700 500]);
% plots
ax1=axes('position',[0.12 0.11 0.8 0.8]);
for k=1:3
bb(k)=plot(2*h0,phixy_optimal{k},'o','linewidth',1.2);hold on;
end
%plot([0 kN*(2*pi)],[0 0],'-k','linewidth',.8);
%plot([0 kN*(2*pi)],[0.65*pi*180/pi 0.65*pi*180/pi],'-k','linewidth',.8);
%hold on; plot([0.44 0.44],[-180,180],'k--');
xlabel('$\ell_{lu}$','interpreter','latex','fontsize',12);
ylabel('Phase $\bar{\hat{\phi}}_{lu}$ (deg)','interpreter','latex','fontsize',12);
set(gca,'ytick',[-180:30:180], 'YScale','linear','XScale','linear');
set(gca,'xminortick','on','yminortick','on',...
    'ticklength',[0.015 0.015],'fontsize',12)
hold on; plot([0.6585*2 0.6585*2],[0 90],'k--');
%ax2.YAxis.MinorTickValues = -180:30:180;
%ylim([-180,180]);
%xlim([0 3]);
legend(bb,'t=26.67','t=53.62','t=80.32','fontsize',16,'location','southeast');

