% This code is to calculate the mixing efficiency, and should be
% run after readmean_h5.m
% clear all;
% User settings....
% Set the run director
% base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/3D_med_5';
addpath(genpath('/glade/u/home/liuchihl/TOOLS'));
NY=721;%1729 % Here, NY should match the value in grid_def.all
N_TH=1; % The number of scalars
Re = 1000; NU=1/Re; % Enter the Reynolds number or viscosity from input.dat
Pr=1;   kappa=NU/Pr; % Prandtl number
Ri_b(1:N_TH)=0.16; % Enter the richardson number for each scalar
fname1 = 'Large_case/d_10/Ri_0.16';fname2 = 'Large_case/d_5/Ri_0.12';
fname3 = 'Large_case/d_3/Ri_0.12';
Lz = 19.86; %vertical length


% Preallocate variables for speed
epsilon=zeros(NY,nk,N_TH);pe_diss=zeros(NY,nk,N_TH);ume=zeros(NY,nk);
thv = zeros(NY,nk,N_TH);dthdy= zeros(NY,nk,N_TH);
gamma1 = []; gamma2 = []; eta = [];EPS=[];
PE_a3=[];PE_b3=[];PE_t3=[];KE_3=[];
for j = 1 % 3 cases
    if j==1
base_dir=['/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/',fname1];
    elseif j==2    
base_dir=['/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/',fname2];
    else
base_dir=['/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/',fname3];
    end
    % Set the filename
filename_tke=[base_dir '/tke.h5'];
filename_mean=[base_dir '/mean.h5'];
filename_buoypdf=[base_dir '/buoypdf.h5'];

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

% Here, read in the statistics calculated in save_stats
% this is from tke.h5, including epsilon; pe_diss is from mean.h5
 
varname=['/epsilon/' timename]; % TKE dissipation rate
epsilon(:,k)=h5read(filename_tke,varname);
varname=['/pe_diss' num2str(n,'%2.2d') '/' timename]; % POTENTIAL ENERGY DISSIPATION
pe_diss(:,k)=h5read(filename_mean,varname);
varname=['/thv' num2str(n,'%2.2d') '/' timename]; % theta*w
thv(:,k)=h5read(filename_mean,varname);
varname=['/dthdy' num2str(n,'%2.2d') '/' timename]; %mean buoyancy gradient 
dthdy(:,k)=h5read(filename_mean,varname);
varname=['/thme' num2str(n,'%2.2d') '/' timename];  % MEAN BUOYANCY
thme(:,k)=h5read(filename_mean,varname);
varname=['/ume/' timename];  % MEAN velocity
ume(:,k)=h5read(filename_mean,varname);

varname=['/shear/' timename];           % MEAN SQUARE SHEAR
shear(:,k)=h5read(filename_mean,varname);

% Read buoyancy pdf files
    if k<10
        binname=['/bins/000' num2str(k)];
        dbinname=['/dbins/000' num2str(k)];
        countname=['/counts/000' num2str(k)];
    elseif (k>=10)&&(k<100)
        binname=['/bins/00' num2str(k)];
        dbinname=['/dbins/00' num2str(k)];
        countname=['/counts/00' num2str(k)];        
    elseif (k>=100)&&(k<1000)
        binname=['/bins/0' num2str(k)];
        dbinname=['/dbins/0' num2str(k)];
        countname=['/counts/0' num2str(k)];  
    else
        binname=['/bins/' num2str(k)];
        dbinname=['/dbins/' num2str(k)];
        countname=['/counts/' num2str(k)];
    end
    
    bins=h5read(filename_buoypdf,binname);
    dbins=h5read(filename_buoypdf,dbinname);
    counts=h5read(filename_buoypdf,countname);
%     time(k)=h5readatt(filename_buoypdf,binname,'Time');
    
    binval=bins+dbins/2;
    
    normalization_factor=trapz(binval,counts);
    buoyancy_pdf(:,k)=counts/normalization_factor;

Z=Lz*cumtrapz(binval,buoyancy_pdf(:,k)); 

% calculate total PE, BPE, APE and KE
KE(k)=.5*mean(((ume(:,k)).^2));           %Volume-averaged kinetic Energ)y
PE_t(k)=-Ri_b*mean(thme(:,k).*gyf');    %Total PE
PE_b(k)=-Ri_b*trapz(Z,(binval.*Z))/Lz;  %Background PE
PE_a(k)=PE_t(k)-PE_b(k);                %Available PE
%
D_pt(k)=Ri_b*(thme(end,k)-thme(1,k))/(Re*Pr*Lz)*time(k); 

end
M(j,:)= mmderiv(time,PE_b-D_pt);
KE_3=[KE_3;KE]; PE_t3=[PE_t3;PE_t]; PE_b3=[PE_b3;PE_b]; PE_a3=[PE_a3;PE_a];
% Brunt-Vaisalla (buoyancy) frequency
% for k=1:nk
%   for j=1:NY
% brunt(j,k)=sqrt((dthdy(j,k))); 
% % Gradient richardson number
% grarich(j,k)=brunt(j,k).^2./shear(j,k); 
% 
%   end
% end
% epsilon_rho = kappa_t.*dthdy.^2;
% epsilon_PE = -RI*thv;
% B = -RI*thv; %buoyancy flux: including the reversible mixing
% clear gamma1 gamma2 gamma1_m gamma2_m
for tt=1:length(time)
% a = find(~isnan(pe_diss(:,tt))&pe_diss(:,tt)<100);
%     pe_diss_m(tt) = mean(pe_diss(a,tt));% avg over depth
b = find(~isnan(epsilon(:,tt)));
    epsilon_m (tt) = mean(epsilon(b,tt));% avg over depth
    epsilon_m(epsilon>1)=nan; %eliminate unreasonble values
% c = find(~isnan(B(:,tt)));
%     B_m(tt) = mean(B(c,tt));
% d = find(~isnan(epsilon_PE(:,tt)));
%     epsilon_PE_m(tt) = mean(epsilon_PE(d,tt));
end
EPS=[EPS;epsilon_m];%TKE dissipation rate
eta = [eta;M(j,:)./(M(j,:)+EPS(j,:))];

% epsilon_PE = RI/(Re*Pr)./pe_diss_m;
% gamma1 = [gamma1;epsilon_PE_m./(epsilon_PE_m+epsilon_m)]; %mixing efficiency:potential energy dissipation
% gamma1 = [gamma1;pe_diss_m./(pe_diss_m+epsilon_m)]; %mixing efficiency:potential energy dissipation
% gamma2 = [gamma2;B_m./epsilon_m];%another definition of mixing efficiency
% clear B_m epsilon_m pe_diss_m
end
%%

close all;
color=[84, 68, 187;213, 21, 21;255, 128, 0]/255;
figure('position',[10 20 500 480]);
ax1 = axes('position',[.17 .71 .8 .26]);
cases=1;
for j=1:cases
    plot(time,KE_3(j,:)-KE_3(j,1),'linewidth',1.5,'color',color(j,:));hold on;
end
set(gca,'fontsize',12,'xticklabel',[]);ylabel('$$\Delta K$$','fontsize',12,'interpreter','latex');
le=legend('d=10','d=5','d=3');set(le,'fontsize',12,'location','southwest');
set(gca,'tickdir','out','xminortick','on','yminortick','on','fontsize',12);

ax2 = axes('position',[.17 .43 .8 .26]);
for j=1:cases
    plot(time,PE_t3(j,:)-PE_t3(j,1),'linewidth',1.5,'color',color(j,:));hold on;
    plot(time,PE_a3(j,:)-PE_a3(j,1),'--','linewidth',1.5,'color',color(j,:));hold on;
end
set(gca,'fontsize',12,'xticklabel',[]);set(gca,'tickdir','out','xminortick','on','yminortick','on','fontsize',12);

ylabel('$$\Delta P_{t}, \Delta P_{a}$$','fontsize',12,'interpreter','latex');

ax3 = axes('position',[.17 .15 .8 .26]);
for j=1:cases
    plot(time,PE_b3(j,:)-PE_b3(j,1)-D_pt,'linewidth',1.5,'color',color(j,:)); hold on;
end
set(gca,'fontsize',12);
set(gca,'tickdir','out','xminortick','on','yminortick','on','fontsize',12,...
    'yticklabel',{'0','0.005','0.01','0.015'});

ylabel('$$\Delta P_{b}-D_{p}t$$','fontsize',12,'interpreter','latex');
xlabel('Time','fontsize',12);
% print -djpeg Energetics.jpg
%%
figure;
set(gcf,'units','centimeters','paperunits','centimeters')
set(gcf,'PaperType','A4');
pp=[0.63 0.9 19.7 10];
ps=[0 0 pp(3)/1.1 pp(4)/1.1];
set(gcf,'paperposition',pp);
set(gcf,'position',ps);
color=[84, 68, 187;213, 21, 21;255, 128, 0]/255;

ax1 = axes('position',[.1 .56 .85 .35]);
for i=1:cases
    
    plot(time,EPS(i,:),'--','linewidth',1.5,'color',color(i,:));
    hold on;
    a(i)=plot(time,M(i,:),'-','linewidth',1.5,'color',color(i,:));
end
xlim([0 time(end)]);
% le=legend([a(1),a(2),a(3)],'d=10','d=5','d=3');
% set(le,'location','northwest');
ylabel('M, \epsilon');
set(gca,'xticklabel',[],'tickdir','out','xminortick','on','yminortick','on','fontsize',12);

ax2 = axes('position',[.1 .18 .85 .35]);
for i=1:cases
plot(time,eta(i,:),'linewidth',1.5,'color',color(i,:));xlim([0 time(end)]);
hold on;
end
ylabel('$$\eta_{i}$$','fontsize',12,'interpreter','latex')
xlabel('Time','fontsize',12);
set(gca,'tickdir','out','xminortick','on','yminortick','on','fontsize',12);
ylim([0 1.5]);
% print -djpeg M_eps_mixginefficenct.jpg
%%
% close all;
% figure; 
% subplot 121;
% for i = 1:3
%     plot(time,gamma1(k,:),'linewidth',2);xlabel('Time','fontsize',12);
% hold on;
% end
% hold on; plot([time(1),time(end)],ones(1,2)*.2,'--k','linewidth',2);
% title('$$\Gamma = P_\epsilon/(P_\epsilon+\epsilon)$$','fontsize',12,'interpreter','latex');
% legend('d_b = 10','d_b = 5','d_b = 3');
% ylim([-1 3])
% subplot 122;
% for i= 1:3
%     plot(time,gamma2(k,:),'linewidth',2);xlabel('Time','fontsize',12);
% hold on;
% end
% hold on; plot([time(1),time(end)],ones(1,2)*.2,'--k','linewidth',2);
% title('$$\Gamma = B/\epsilon$$','fontsize',12,'interpreter','latex');
% set(gcf,'color','w');


%%
% pe_diss(pe_diss>100) = nan;
% % Rf = pe_diss./(pe_diss+epsilon);
% epsilon_PE = shear.^2.*kappa_t*RI;
% Rf = epsilon_PE./(epsilon_PE+epsilon);
% % Rf(Rf>1|Rf<-1)=nan;
% figure; 
% figset(2,1);
% axes('position',[.1 .55 .4 .4]);
% pcolor(time,gyf,log10(epsilon));shading interp;
% caxis([-8 -3]);
% colorbar;
% title('log\epsilon','fontsize',12);
% axes('position',[.52 .55 .4 .4]);
% pcolor(time,gyf,epsilon_PE);shading interp;
% colorbar;
% title('\epsilon_P_E','fontsize',12);
% axes('position',[.1 .1 .4 .4]);
% pcolor(time,gyf,(Rf));shading interp;
% colorbar;
% title('R_f','fontsize',12);


