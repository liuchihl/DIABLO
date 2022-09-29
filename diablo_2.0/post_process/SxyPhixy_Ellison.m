% This script calculates the cross-spectra between upper and lower disturbances without smoothing the speectrum via the Ellison scale instead of the vertical velocity perturbation (v')
% Here we can either use 2D slices or 3D flow field
% see crossspec.m for squared coherence (smoothed), see 
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

for k=1:2

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

varname=['/Timestep/TH1'];

if exist([filename,timename])
TH1(:,:,:)=h5read([filename,timename],varname);

info = h5info([filename,timename]);
%time
tt(k) = info.Groups.Attributes.Value;

% Ellison Scale LE: for consistency with the textbook (Smyth and Carpenter) Figure 3.22, instead using the vertical velocity perturbation, the vertical displacement is perhaps good to try
B = mean(TH1,[1,3]); 		  % horizontal averaged buoyancy
By = gradient(B,mean(diff(y)));   % vertical buoayncy gradient 
b = TH1-B;			  % root-mean-square of buoyancy perturbation
LE = b./repmat(By,[NX,1,NZ]);	  % Ellison scale


% find the upper and lower layer that need to be crosscorrelated
[~,up]=min(abs(y-0.6585));
[~,low]=min(abs(y+0.6585));
LEup = squeeze(LE(:,up,:)); LElow = squeeze(LE(:,low,:));
N=length(LEup);
% compute the cross correlation with the written function: fft_cpsd
for j=1:nz
	[Sxx,Syy,Sxy,kj] = fft_cpsd(LElow(:,j),LEup(:,j),dx,'rec');

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
Sxyk{k} = temp1./nz;
phixyk{k} = temp2./nz;
end
end

%close all;
hh=figure('position',[50 50 700 800]);
% plots
ax1=axes('position',[0.12 0.52 0.8 0.42]);
for kk=1:3
aa(kk)=semilogy(kj.*(2*pi),real(Sxyk{kk}),'-','linewidth',1.2); hold on;
end
xlabel('$k^\ast$','interpreter','latex','fontsize',12);
ylabel('CPSD','interpreter','latex','fontsize',12);
hold on; 
hold on; plot([0.44 0.44],[0,1],'k--');
set(gca, 'YScale','linear','XScale','linear','xminortick','on','yminortick','on',...
    'ticklength',[0.015 0.015],'fontsize',12,'xticklabel',[])
xlim([0 3]) ;ylim([0 0.04]);
legend(aa,'t=26.67','t=53.62','t=80.32','fontsize',16);

ax2=axes('position',[0.12 0.09 0.8 0.42]);
for kk=1:3
bb=plot(kj*2*pi,phixyk{kk},'-o','linewidth',1.2);hold on;
end
plot([0 kN*(2*pi)],[0 0],'-k','linewidth',.8);
plot([0 kN*(2*pi)],[0.65*pi*180/pi 0.65*pi*180/pi],'-k','linewidth',.8);
hold on; plot([0.44 0.44],[-180,180],'k--');
xlabel('$k^\ast$','interpreter','latex','fontsize',12);
ylabel('Phase $\bar{\hat{\phi}}_{lu}(k^\ast)$ (deg)','interpreter','latex','fontsize',12);
set(gca,'ytick',[-180:120:180], 'YScale','linear','XScale','linear');
set(gca,'xminortick','on','yminortick','on',...
    'ticklength',[0.015 0.015],'fontsize',12)
ax2.YAxis.MinorTickValues = -180:30:180;
ylim([-180,180]);
xlim([0 3]);


