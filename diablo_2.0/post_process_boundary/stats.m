 % This script computes stats for butterfly effects and save data into .mat

% Run after readmean.m
%clear;%clc; 
addpath(genpath('/glade/u/home/liuchihl/TOOLS'));
base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_2.5/Ri_0.12_0.1_';

% input coefficients

NY=361;LY=20;
N_TH=1;
Re = 1000; NU=1/Re; % Enter the Reynolds number or viscosity from input.dat
Pr=1;   kappa=NU/Pr; % Prandtl number
RI(1:N_TH)=0.12; % Enter the richardson number for each scalar
num=1;		% number of ensembles
nk = zeros(1,num);
for i=1:num
% Set the filename
filename=[base_dir num2str(i) '_rightBC/mean.h5'];
filename_tke=[base_dir num2str(i) '_rightBC/tke.h5'];
filename_buoypdf=[base_dir num2str(i) '_rightBC/buoypdf.h5'];

file_info=h5info(filename);
att_info=file_info.Groups.Attributes;
nk(i)=att_info.Value
% Preallocate variables for speed
time=zeros(1,nk(i));
ume{i}=zeros(NY,nk(i)); vme=zeros(NY,nk(i)); wme=zeros(NY,nk(i));
urms=zeros(NY,nk(i)); vrms=zeros(NY,nk(i)); wrms=zeros(NY,nk(i));
uv=zeros(NY,nk(i)); uw=zeros(NY,nk(i)); wv=zeros(NY,nk(i));
dudy=zeros(NY,nk(i)); dvdy=zeros(NY,nk(i)); dwdy=zeros(NY,nk(i));
cp=zeros(NY,nk(i)); shear=zeros(NY,nk(i));
omega_x=zeros(NY,nk(i)); omega_y=zeros(NY,nk(i)); omega_z=zeros(NY,nk(i));
thme{i}=zeros(NY,nk(i),N_TH);dthdy=zeros(NY,nk(i),N_TH); thrms=zeros(NY,nk(i),N_TH);
Time{i}=zeros(1,nk(i));

thv=zeros(NY,nk(i),N_TH); pe_diss=zeros(NY,nk(i),N_TH); k2d{i}=zeros(NY,nk(i)); k3d{i}=zeros(NY,nk(i));
epsilon{i}=zeros(NY,nk(i)); Dpt{i}=zeros(1,nk(i));
B_flux{i}=zeros(NY,nk(i)); shear_prod{i}=zeros(NY,nk(i));
EPS{i}=zeros(NY,nk(i)); M{i}=zeros(NY,nk(i)); eta{i}=zeros(NY,nk(i));
Mixingrate=zeros(1,nk(i));
PE_b=zeros(1,nk(i)); PE_t=zeros(1,nk(i)); PE_a=zeros(1,nk(i)); KE_t=zeros(1,nk(i));
KE_3{i}=zeros(1,nk(i)); D_pt=zeros(1,nk(i)); PE_t3{i}=zeros(1,nk(i)); PE_b3{i}=zeros(1,nk(i));
PE_a3{i}=zeros(1,nk(i)); epsilon_m=zeros(1,nk(i)); tke{i}=zeros(NY,nk(i));
N2{i}=zeros(NY,nk(i));

for k=1:nk(i)



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
    time(k)=h5read(filename,varname);
    varname=['/gyf/' timename];             % Y-COORDINATE
    gyf(:)=h5read(filename,varname);

    varname=['/ume/' timename];             % MEAN VELOCITIES
    ume{i}(:,k)=h5read(filename,varname);
    varname=['/vme/' timename];
    vme(:,k)=h5read(filename,varname);
    varname=['/wme/' timename];
    wme(:,k)=h5read(filename,varname);
    varname=['/urms/' timename];            % RMS VELOCITIES
    urms(:,k)=h5read(filename,varname);
    varname=['/vrms/' timename];
    vrms(:,k)=h5read(filename,varname);
    varname=['/wrms/' timename];
    wrms(:,k)=h5read(filename,varname);
    varname=['/uv/' timename];              % REYNOLDS STRESSES
    uv(:,k)=h5read(filename,varname);
    varname=['/uw/' timename];
    uw(:,k)=h5read(filename,varname);
    varname=['/wv/' timename];
    wv(:,k)=h5read(filename,varname);
    varname=['/dudy/' timename];            % MEAN VELOCITY GRADIENTS

    dwdy(:,k)=h5read(filename,varname);
    varname=['/shear/' timename];           % MEAN SQUARE SHEAR
    shear(:,k)=h5read(filename,varname);
    varname=['/omega_x/' timename];         % RMS VORTICITIES
    omega_x(:,k)=h5read(filename,varname);
    varname=['/omega_y/' timename];
    omega_y(:,k)=h5read(filename,varname);
    varname=['/omega_z/' timename];
    omega_z(:,k)=h5read(filename,varname);
    varname=['/k2d/' timename];
    k2d{i}(:,k)=h5read(filename,varname);
    varname=['/k3d/' timename];
    k3d{i}(:,k)=h5read(filename,varname);
    varname=['/epsilon/' timename]; % TKE dissipation rate
    epsilon{i}(:,k)=h5read(filename_tke,varname);
    for n=1:N_TH
        varname=['/thme' num2str(n,'%2.2d') '/' timename];  % MEAN BUOYANCY
        thme{i}(:,k,n)=h5read(filename,varname);
        varname=['/dthdy' num2str(n,'%2.2d') '/' timename]; % MEAN BUOYANCY GRADIENTS
        dthdy(:,k,n)=h5read(filename,varname);
        varname=['/thrms' num2str(n,'%2.2d') '/' timename]; % RMS BUOYANCY
        thrms(:,k)=h5read(filename,varname);
        varname=['/thv' num2str(n,'%2.2d') '/' timename]; % theta*w
        thv(:,k)=h5read(filename,varname);
        varname=['/pe_diss' num2str(n,'%2.2d') '/' timename]; % POTENTIAL ENERGY DISSIPATION
        pe_diss(:,k)=h5read(filename,varname);
    end
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
    binval=bins+dbins/2;
    normalization_factor=trapz(binval,counts);
    buoyancy_pdf(:,k)=counts/normalization_factor;
    Y=LY*cumtrapz(binval,buoyancy_pdf(:,k)); 

% calculate total PE, BPE, APE and KE


%Volume-averaged total kinetic Energy
KE_t(k)=.5*mean((urms(:,k).^2+vrms(:,k).^2+...
    wrms(:,k).^2)+(ume{i}(:,k).^2+...
    vme(:,k).^2+wme(:,k).^2));                 %Total KE

%Volume-averaged Potential Energy
PE_t(k)=-RI*mean( thme{i}(:,k)' .*gyf);  %Total PE
PE_b(k)=-RI*trapz(Y,binval.*Y)/LY;    %Background PE
PE_a(k)=PE_t(k)-PE_b(k);              %Available PE
%Volume-averaged of laminar diffusion of mean flow
D_pt(k)=RI*(thme{i}(end,k)-thme{i}(1,k))/(Re*Pr*LY)*time(k); 

end    %(end of time)
% Combine all cases
KE_3{i}=KE_t; PE_t3{i}=PE_t; PE_b3{i}=PE_b;
PE_a3{i}=PE_a;
Dpt{i}=D_pt;
Mixingrate= mmderiv(time,PE_b3{i}-Dpt{i});

% total potential energy dissipation
epsilon_P{i} = pe_diss.*RI./(Re*Pr*dthdy);     
% perturbation potential energy dissipation
epsilon_Pprime{i}=(pe_diss-dthdy.^2)*RI./(Re*Pr*dthdy) ;
epsilon_Pprime{i}(epsilon_Pprime{i}<-20)=0;     % eliminate 1 peak that contaminate our results
B_flux{i} = RI.*mean(thv); 			% buoyancy flux: including the reversible mixing
shear_prod{i} = -mean(uv.*dudy);                % Shear production uw*dU/dz
dvdy = mmderiv(gyf,vme);

%TKE dissipation rate (epsilon_m) + mean flow dissipation rate (epsilon_mf_m)
clear epsilon_m
for tt=1:nk(i)
    epsilon_m(tt) = nanmean(epsilon{i}(3:end,tt));% avg over depth
end
% epsilon_mf_m = 1/Re*(mean(shear));
epsilon_mf_m = 1/Re*(mean(dudy(3:end,:).^2+dvdy(3:end,:).^2+dwdy(3:end,:).^2));
EPS{i}=epsilon_m+epsilon_mf_m;           % total epsilon
M{i}=Mixingrate;
eta{i} = M{i}./(M{i}+EPS{i});            %mixing efficiency eta
EPS{i}(EPS{i}>1)=0;
% time
    Time{i} = time;
% calculate tke
    tke{i}=0.5*(urms.^2.+vrms.^2.+wrms.^2.);
% Gradient Richardson number N2/S2
    Rig{i} = (RI*dthdy)./shear;
% Buoyancy Re ( 1/Ri_b* <dui/dxj*dui/dxj> ): use the total velocity
    Re_b{i} = 1./(RI*dthdy(3:end,:)).*( epsilon{i}(3:end,:)/NU+shear(3:end,:) );
% Buoyancy frequency square
    N2{i} = RI*dthdy;




% terms of integrals related to mixing efficiency

% clear EPS_in M_in Rf Gamma_c 
j=3; 
for i=1:num
    EPSp_in{i}=zeros(1,nk(i));EPS_in{i}=zeros(1,nk(i)); M_in{i}=zeros(1,nk(i));
    eta_c{i}=zeros(1,nk(i));Rf_star{i}=zeros(1,nk(i)); 
%    Time{i} = Time{i}(1:nk(i));
for ff=j+1:nk(i)
    M_in{i}(ff) = trapz(Time{i}(j:ff),M{i}(j:ff));
    EPS_in{i}(ff) = trapz(Time{i}(j:ff),EPS{i}(j:ff)); 
    EPSp_in{i}(ff) = trapz(Time{i}(j:ff),mean(epsilon_Pprime{i}(3:end-1,j:ff)));
end
% definitions of cumulative mixing efficiency
    eta_c{i} = M_in{i}./(M_in{i}+EPS_in{i});         % Cumulative Flux coeeficient
    Gamma_c{i} = M_in{i}./(EPS_in{i}); 
    Rf_star{i} = EPSp_in{i}./(EPS_in{i}+EPSp_in{i});         % Cumulative irreversible flux Richardson number
%Dp_in(ff) = trapz(Time{j}(i:ff),Dpt{j}(i:ff)./Time{j}(i:ff));
end

end
% save data that is useful for Ri comparison

save CumValue_12_25_sta_rightBC M_in EPS_in EPSp_in eta_c Rf_star Time
save stats_12_25_sta_rightBC ume thme k2d k3d KE_3 PE_b3 PE_t3 PE_a3 tke B_flux shear_prod EPS M Dpt epsilon_Pprime epsilon eta Rig Re_b N2 Time 



%% save mean flow



%clear;%clc; 
%addpath(genpath('/glade/u/home/liuchihl/TOOLS'));
%base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_10/butterfly/Ri_0.12_0.05_';

% input coefficients

%NY=361;LY=20;
%N_TH=1;
%Re = 1000; NU=1/Re; % Enter the Reynolds number or viscosity from input.dat
%Pr=1;   kappa=NU/Pr; % Prandtl number
%RI(1:N_TH)=0.12; % Enter the richardson number for each scalar
%num=15; 	 % number of cases
%nk=501;
%for i=1:num
% Set the filename
%filename=[base_dir num2str(i) '/mean.h5'];

% Preallocate variables for speed
%time=zeros(1,ni);
%ume=zeros(NY,ni); vme=zeros(NY,ni); wme=zeros(NY,ni);
%thme=zeros(NY,ni);

%for k=1:nk 
%
%    if (k<10)
%        timename=['000' int2str(k)];
%    elseif (k<100)
%        timename=['00' int2str(k)];
%    elseif (k<1000)
%        timename=['0' int2str(k)];
%    else
%        timename=[int2str(k)];
%    end
%    varname=['/time/' timename];            % TIME
%    time(k)=h5read(filename,varname);
%    varname=['/gyf/' timename];             % Y-COORDINATE
%    gyf(:)=h5read(filename,varname);

%    varname=['/ume/' timename];             % MEAN VELOCITIES
%    ume(:,k)=h5read(filename,varname);
    
%  n=1;
%    varname=['/thme' num2str(n,'%2.2d') '/' timename];  % MEAN BUOYANCY
%    thme(:,k,n)=h5read(filename,varname);
        

%end
%i
%end
%Z=gyf;
%save meanprofile12 time ume thme Z
