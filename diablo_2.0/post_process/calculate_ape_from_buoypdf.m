addpath(genpath('/glade/u/home/liuchihl/TOOLS'));

%gyf=h5read('../grid.h5','/grids/y');
%nt=1001;
%note the coordinate here is a standard Cartesian, i.e.,
%(x,y,z) = (alongstream,cross-stream, vertical)
%Shouldn't be confused with the model definition of Y and Z..
z=(gyf(1:end-1)+gyf(2:end))/2;
nx=384;ny=128;
x=linspace(0,41.88,nx);
y=linspace(0,7.0,ny);
Ri_b=0.15;Re=1000;Pr=1;
Lz=20; % vertical domain

for ts=1:nk

figure(2)
subplot(1,2,1)
plot(binval,buoyancy_pdf(:,ts),'-')
drawnow
ylim([0 1.5]);%xlim([-1.5 1.5]);
Z=Lz*cumtrapz(binval,buoyancy_pdf(:,ts)); 
% 
figure(2)
subplot(1,2,2)
plot(binval,Z)
drawnow
title(sprintf('Time=%0.3f',time(ts)));
% calculate total PE, BPE, and APE
PE_b(ts)=-Ri_b*trapz(Z,(binval.*Z))/Lz;
D_pt(ts)=-Ri_b*(-1-1)/(Re*Pr*Lz)*time(ts); % say for now delta(t) is approximately 5e-2

end
M= mmderiv(time,PE_b-D_pt);
figure; hold on;plot(time,M);
legend('d=10','d=5','d=3');set(gca,'fontsize',12);
xlabel('Time');ylabel('$$M(dP_{b}/dt-D_{p})$$','fontname','times','fontangle','italic','interpreter','latex')
figure;
plot(time,PE_b-PE_b(1)-D_pt,'linewidth',2)
legend('PE_b-D_pt-PE_b(0)')

%save ape_Re1000Pr1Ri012_sk01 PE_b D_pt buoyancy_pdf binval x y z time