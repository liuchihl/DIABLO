 % This script computes stats and save data into .mat

% Run after readmean.m
%clear;%clc; 
addpath(genpath('/glade/u/home/liuchihl/TOOLS'));
for jj=[1:6]
if jj==1
base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/doubleshearlayer/Ri016/R1/3D/D_0.5_';
elseif jj==2
base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/doubleshearlayer/Ri016/R1/3D/D_1_'
elseif jj==3
base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/doubleshearlayer/Ri016/R1/3D/D_1.5_'
elseif jj==4
base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/doubleshearlayer/Ri016/R1/3D/D_2_'
elseif jj==5
base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/doubleshearlayer/Ri016/R1/3D/D_3_';
else jj==6
base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/doubleshearlayer/Ri016/R1/3D/D_0_';
end
% input coefficients

NY=613;LY=30;
N_TH=1;
Re = 1000; NU=1/Re; % Enter the Reynolds number or viscosity from input.dat
Pr=1;   kappa=NU/Pr; % Prandtl number
RI(1:N_TH)=0.16; % Enter the richardson number for each scalar

cs = [1:3];
num=length(cs);		% number of ensembles
nk = zeros(1,num);
for i=1:num
% Set the filename
filename=[base_dir num2str(cs(i)) '/mean.h5'];
filename_tke=[base_dir num2str(cs(i)) '/tke.h5'];
filename_buoypdf=[base_dir num2str(cs(i)) '/buoypdf.h5'];

file_info=h5info(filename);
att_info=file_info.Groups.Attributes;
nk(i)=att_info.Value
clear time
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
end

% eliminate time discontinuity result from multiple runs
k1 = 1:nk(i);
while any(diff(time)<0.1)
ind1 = find(diff(time)<0.1);
k1(ind1(1)+1) = [];
time(ind1(1)+1) = [];
end

% Preallocate variables for speed
nk(i) = length(time);  % redefine nk since we eliminate some discontinuous data
%time=zeros(1,nk(i));
ume{i}=zeros(NY,nk(i)); vme=zeros(NY,nk(i)); wme=zeros(NY,nk(i));
urms=zeros(NY,nk(i)); vrms=zeros(NY,nk(i)); wrms=zeros(NY,nk(i));
uv=zeros(NY,nk(i)); uw=zeros(NY,nk(i)); wv=zeros(NY,nk(i));
dudy=zeros(NY,nk(i)); dvdy=zeros(NY,nk(i)); dwdy=zeros(NY,nk(i));
cp=zeros(NY,nk(i)); shear{i}=zeros(NY,nk(i));
omega_x=zeros(NY,nk(i)); omega_y=zeros(NY,nk(i)); omega_z=zeros(NY,nk(i));
thme{i}=zeros(NY,nk(i),N_TH);dthdy=zeros(NY,nk(i),N_TH); thrms=zeros(NY,nk(i),N_TH);
Time{i}=zeros(1,nk(i));

thv=zeros(NY,nk(i),N_TH); pe_diss=zeros(NY,nk(i),N_TH); k2d{i}=zeros(NY,nk(i)); k3d{i}=zeros(NY,nk(i));
epsilon{i}=zeros(NY,nk(i)); Dpt{i}=zeros(NY,nk(i));
B_flux{i}=zeros(NY,nk(i)); shear_prod{i}=zeros(NY,nk(i));
EPS{i}=zeros(NY,nk(i)); M{i}=zeros(NY,nk(i)); eta{i}=zeros(NY,nk(i));
Mixingrate=zeros(NY,nk(i));
PE_b=zeros(NY,nk(i)); PE_t=zeros(NY,nk(i)); PE_a=zeros(NY,nk(i)); KE_t=zeros(NY,nk(i));
KE_3{i}=zeros(NY,nk(i)); D_pt=zeros(NY,nk(i)); PE_t3{i}=zeros(NY,nk(i)); PE_b3{i}=zeros(NY,nk(i));
PE_a3{i}=zeros(NY,nk(i)); epsilon_m=zeros(NY,nk(i)); tke{i}=zeros(NY,nk(i)); dbstdy{i}=zeros(1,nk(i));
N2{i}=zeros(NY,nk(i));

for k=1:length(k1)-1
k
    if (k1(k)<10)
        timename=['000' int2str(k1(k))];
    elseif (k1(k)<100)
        timename=['00' int2str(k1(k))];
    elseif (k1(k)<1000)
        timename=['0' int2str(k1(k))];
    else
        timename=[int2str(k1(k))];
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
    dudy(:,k)=h5read(filename,varname);
    varname=['/dwdy/' timename];            % MEAN VELOCITY GRADIENTS
    dwdy(:,k)=h5read(filename,varname);
    varname=['/shear/' timename];           % MEAN SQUARE SHEAR
    shear{i}(:,k)=h5read(filename,varname);
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
    if k1(k)<10
        binname=['/bins/000' num2str(k1(k))];
        dbinname=['/dbins/000' num2str(k1(k))];
        countname=['/counts/000' num2str(k1(k))];
    elseif (k1(k)>=10)&&(k1(k)<100)
        binname=['/bins/00' num2str(k1(k))];
        dbinname=['/dbins/00' num2str(k1(k))];
        countname=['/counts/00' num2str(k1(k))];
    elseif (k1(k)>=100)&&(k1(k)<1000)
        binname=['/bins/0' num2str(k1(k))];
        dbinname=['/dbins/0' num2str(k1(k))];
        countname=['/counts/0' num2str(k1(k))];
    else
        binname=['/bins/' num2str(k1(k))];
        dbinname=['/dbins/' num2str(k1(k))];
        countname=['/counts/' num2str(k1(k))];
    end
    bins=h5read(filename_buoypdf,binname);
    dbins=h5read(filename_buoypdf,dbinname);
    counts=h5read(filename_buoypdf,countname);
    binval=bins+dbins/2;
    normalization_factor=trapz(binval,counts);
    buoyancy_pdf(:,k)=counts/normalization_factor;
    Y=LY*cumtrapz(binval,buoyancy_pdf(:,k)); 

% calculate total PE, BPE, APE and KE


%total kinetic Energy
KE_t(:,k)=.5*((urms(:,k).^2+vrms(:,k).^2+...
    wrms(:,k).^2)+(ume{i}(:,k).^2+...
    vme(:,k).^2+wme(:,k).^2));                 %Total KE

%Potential Energy
PE_t(:,k)=-RI*( thme{i}(:,k)' .*gyf);  %Total PE
PE_b(:,k)=-RI*trapz(Y,binval.*Y)/LY;    %Background PE
PE_a(:,k)=PE_t(:,k)-PE_b(:,k);              %Available PE

% db*/dz (b* is the sorted buoyancy profile)
temp = mmderiv(Y(binval<1 & Y>0),binval(binval<1 & Y>0));
%lm = max(temp(100:end));   % find the maximum value away from the boundary, and ignore values that are larger than that
dbstdy{i}(k)=mean(temp(abs(temp)<2));  % elimnate large derivatives near boundary, 2 is arbitrarily chosen
%Volume-averaged of laminar diffusion of mean flow
thme_y = mmderiv(gyf,thme{i}(:,k));
thme_yy = mmderiv(gyf,thme_y);
D_pt(:,k)=-RI*gyf'.*thme_yy./(Re*Pr)*time(k);
%D_pt(k)=RI*(thme{i}(end,k)-thme{i}(1,k))/(Re*Pr*LY)*time(k); 

end    %(end of time)
% Combine all cases
KE_3{i}=KE_t; PE_t3{i}=PE_t; PE_b3{i}=PE_b;
PE_a3{i}=PE_a;
Dpt{i}=D_pt;
Mixingrate= mmderiv(time,PE_b3{i}'-Dpt{i}')';

% total potential energy dissipation
epsilon_P{i} = pe_diss.*RI./(Re*Pr*dthdy);     
% perturbation potential energy dissipation
epsilon_Pprime{i}=(pe_diss-dthdy.^2)*RI./(Re*Pr*dthdy) ;
epsilon_Pprime{i}(epsilon_Pprime{i}<-20)=0;     % eliminate 1 peak that contaminate our results
B_flux{i} = RI.*(thv); 			% buoyancy flux: including the reversible mixing
shear_prod{i} = -(uv.*dudy);                % Shear production uw*dU/dz
dvdy = mmderiv(gyf,vme);

%TKE dissipation rate (epsilon_m) + mean flow dissipation rate (epsilon_mf_m)
clear epsilon_m
for tt=1:length(k1)
    epsilon_m(tt) = nanmean(epsilon{i}(2:end,tt));% avg over depth
end
% epsilon_mf_m = 1/Re*(mean(shear));
epsilon_mf_m = 1/Re*(mean(dudy(2:end,:).^2+dvdy(2:end,:).^2+dwdy(2:end,:).^2));
EPS{i}=epsilon_m+epsilon_mf_m;           % total epsilon
epsilon_t{i} = epsilon{i}+1/Re*(dudy.^2+dvdy.^2+dwdy.^2);
%EPS{i}=epsilon_m;           % epsilon'
M{i}=Mixingrate;
eta{i} = M{i}./(M{i}+epsilon{i});            %mixing efficiency eta
eta_avg{i} = mean(M{i})./(mean(M{i})+(EPS{i}));
EPS{i}(EPS{i}>1)=0;
% time
    Time{i} = time;
% calculate tke
    tke{i}=0.5*(urms.^2.+vrms.^2.+wrms.^2.);
% Gradient Richardson number N2/S2
%    Rig{i} = (RI*dthdy)./shear{i};
    Rig{i} = RI*dthdy./(dudy.^2+dwdy.^2);
% Buoyancy Re ( 1/Ri_b* <dui/dxj*dui/dxj> ): use the total velocity
%    Re_b{i} = 1./mean(RI*dthdy(2:end,:)).*( (epsilon_mf_m+epsilon_m)/NU);
     Re_b{i} = 1./(RI*dthdy(1:end,:)).*epsilon_t{i};
% Buoyancy frequency square
    N2{i} = RI*dthdy;
% shear square
    S2{i} = dudy.^2+dwdy.^2;



% terms of integrals related to mixing efficiency

% clear EPS_in M_in Rf Gamma_c 
j=3; 
% for i=1:num
    EPSp_in{i}=zeros(1,nk(i));EPS_in{i}=zeros(1,nk(i)); M_in{i}=zeros(1,nk(i));
    eta_c{i}=zeros(1,nk(i));Rf_star{i}=zeros(1,nk(i)); 
%    Time{i} = Time{i}(1:nk(i));
for ff=j+1:length(k1)
    M_in{i}(ff) = trapz(Time{i}(j:ff),mean(M{i}(:,j:ff)));
    EPS_in{i}(ff) = trapz(Time{i}(j:ff),EPS{i}(j:ff)); 
    EPSp_in{i}(ff) = trapz(Time{i}(j:ff),nanmean(epsilon_Pprime{i}(3:end-1,j:ff)));
end
% definitions of cumulative mixing efficiency
    eta_c{i} = M_in{i}./(M_in{i}+EPS_in{i});         % Cumulative Flux coeeficient
    Gamma_c{i} = M_in{i}./(EPS_in{i}); 
    Rf_star{i} = EPSp_in{i}./(EPS_in{i}+EPSp_in{i});         % Cumulative irreversible flux Richardson number
%Dp_in(ff) = trapz(Time{j}(i:ff),Dpt{j}(i:ff)./Time{j}(i:ff));
end

% end
% save data that is useful for Ri comparison
if jj==1
save CumValue_16_0.5D.mat M_in EPS_in EPSp_in eta_c Rf_star Time
save stats_16_0.5D.mat ume thme k2d k3d KE_3 PE_b3 PE_t3 PE_a3 tke B_flux shear_prod EPS epsilon_P M Dpt epsilon_Pprime epsilon epsilon_t eta Rig Re_b dbstdy N2 S2 Time eta_avg shear gyf
elseif jj==2
save CumValue_16_D1.mat M_in EPS_in EPSp_in eta_c Rf_star Time
save stats_16_D1.mat ume thme k2d k3d KE_3 PE_b3 PE_t3 PE_a3 tke B_flux shear_prod EPS epsilon_P M Dpt epsilon_Pprime epsilon epsilon_t eta Rig Re_b dbstdy N2 S2 Time eta_avg shear gyf
elseif jj==3
save CumValue_16_D1.5.mat M_in EPS_in EPSp_in eta_c Rf_star Time
save stats_16_D1.5.mat ume thme k2d k3d KE_3 PE_b3 PE_t3 PE_a3 tke B_flux shear_prod EPS epsilon_P M Dpt epsilon_Pprime epsilon epsilon_t eta Rig Re_b dbstdy N2 S2 Time eta_avg shear gyf
elseif jj==4
save CumValue_16_D2.mat M_in EPS_in EPSp_in eta_c Rf_star Time
save stats_16_D2.mat ume thme k2d k3d KE_3 PE_b3 PE_t3 PE_a3 tke B_flux shear_prod EPS M Dpt epsilon_P epsilon_Pprime epsilon epsilon_t eta Rig Re_b dbstdy N2 S2 Time eta_avg shear gyf
elseif jj==5
save CumValue_16_D3.mat M_in EPS_in EPSp_in eta_c Rf_star Time
save stats_16_D3.mat ume thme k2d k3d KE_3 PE_b3 PE_t3 PE_a3 tke B_flux shear_prod EPS M Dpt epsilon_P epsilon_Pprime epsilon epsilon_t eta Rig Re_b dbstdy N2 S2 Time eta_avg shear gyf
else jj==6
save CumValue_16_D0.mat M_in EPS_in EPSp_in eta_c Rf_star Time
save stats_16_D0.mat ume thme k2d k3d KE_3 PE_b3 PE_t3 PE_a3 tke B_flux shear_prod EPS M Dpt epsilon_P epsilon_Pprime epsilon epsilon_t eta Rig Re_b dbstdy N2 S2 Time eta_avg shear gyf
end

end
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
