% This script illustrates how to load in 3D model data and calculate 
% some basic diagnostics
% Run after readmean.m


fname1 = 'Ri_0.16_0.5_20_small';
%fname2 = '3D_med_7_Lz6';
%fname3 = '3D_med_7_Lz7';
% fname2 = '3D_med_Ly5';
% fname3 = '3D_med_Ly6';fname4 = '3D_med_Ly7';
% fname5 = '3D_med_Ly8';
LX=27.73;
LZ=6.93;
%inv = 4000;
timestep = [2000:2000:8000];
tt = zeros(3,length(timestep));
tke_3D_m = zeros(1,length(timestep));
for i = 1 %5 different cases
if     i==1    
    filename=['/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_10/' fname1 '/output/'];
% elseif i==2
%     filename=['/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/' fname2 '/output/'];
% else 
%     filename=['/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/' fname3 '/output/'];
% elseif i==4
%     filename=['/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/' fname4 '/output/'];
% else
%     filename=['/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/' fname5 '/output/'];
end

for k=timestep

  if (k<10)
    timename=['out0000' int2str(k) '.h5'];
  elseif (k<100)
    timename=['out000' int2str(k) '.h5'];
  elseif (k<1000)
    timename=['out00' int2str(k) '.h5'];
  elseif (k<10000)
    timename=['out0' int2str(k) '.h5'];
  else
    timename=['out' int2str(k) '.h5'];
  end
varname1=['/Timestep/U'];
varname2=['/Timestep/V'];
varname3=['/Timestep/W'];
% varname4=['/Timestep/TH1'];
if exist([filename,timename])
U=h5read([filename,timename],varname1);
V=h5read([filename,timename],varname2);
W=h5read([filename,timename],varname3);
% TH1=h5read([filename,timename],varname4);
info = h5info([filename,timename]);
%time
tt(i,k/inv) = info.Groups.Attributes.Value;
NX=size(U,1); NZ=size(U,3); NY=size(U,2);
x=linspace(0,LX,NX); z=linspace(0,LZ,NZ);

drhodz1=0.0;
drhodz2=0.0;
% Add the background buoyancy gradient
% for j=1:NZ
%   TH1(:,:,j)=TH1(:,:,j)+drhodz1*z(j);
% %  TH2(:,:,k)=TH2(:,:,k)+drhodz2*z(k);
% end

% calculate tke 3D
U_mean = mean(U,3); U_rms = sqrt((U-U_mean).^2);
V_mean = mean(V,3); V_rms = sqrt((V-V_mean).^2);
W_mean = mean(W,3); W_rms = sqrt((W-W_mean).^2);
tke_3D = .5*(U_rms.^2+V_rms.^2+W_rms.^2);
tke_3D_m(i,k/inv) = mean(mean(mean(tke_3D)));
end

end
end
%%
 
%plot
figure; 
for i=1
plot(tt(i,:),tke_3D_m(i,:),'o-');
hold on;
end
% le = legend('Ly = 4','Ly = 5','Ly = 6','Ly = 7','Ly = 8');
% le = legend('Ly = 4','Ly = 6','Ly = 7');
% set(le,'fontsize',12);
xlabel('time','fontsize',12); ylabel('K_3_D','fontsize',12);
% % Calculate the x-average velocity
% ume=squeeze(mean(U,1));
% vme=squeeze(mean(V,1));
% wme=squeeze(mean(W,1));
% thme1=squeeze(mean(TH1,1));
% %thme2=squeeze(mean(TH2,1));
% 
% % Calculate correlation terms
% for k=1:NZ
%   uw_mean(:,k)=mean((U(:,:,k)-U_BT(k)).*(W(:,:,k)-W_BT(k)),1);
%   uu_mean(:,k)=mean((U(:,:,k)-U_BT(k)).*(U(:,:,k)-U_BT(k)),1);
%   ww_mean(:,k)=mean((W(:,:,k)-W_BT(k)).*(W(:,:,k)-W_BT(k)),1);
%   uw_BT(:,k)=trapz(gyf,uw_mean(:,k))/(gyf(end)-gyf(1));
%   uu_BT(:,k)=trapz(gyf,uu_mean(:,k))/(gyf(end)-gyf(1));
%   ww_BT(:,k)=trapz(gyf,ww_mean(:,k))/(gyf(end)-gyf(1));
% end
% 
% % calculate the mean shear
% dudy=zeros(size(ume));
% dwdy=zeros(size(wme));
% for j=2:NY-1
%     dudy(j,:)=(ume(j+1,:)-ume(j-1,:))/(gyf(j+1)-gyf(j-1));
%     dwdy(j,:)=(wme(j+1,:)-wme(j-1,:))/(gyf(j+1)-gyf(j-1));
% end
% 
% % Calculate the local shear
% dUdy=zeros(size(U));
% for j=2:NY-1
% dUdy(:,j,:)=(U(:,j+1,:)-U(:,j-1,:))/(gyf(j+1)-gyf(j-1));
% end


