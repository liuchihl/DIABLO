% This script is part of the boundary proximity project: it calculates the  phase difference between the KH mode and subharmonic mode at the middle of the shear layer
% see crossspec.m for squared coherence (smoothed)
% In this section, data is being extracted and variables are calculated and
% saved
clear;%clc;

% Define file names and directory
%fname1 = 'Ri_0.16_0.5_20_small/';
 base_dir='/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/d_10/boundary/Ri_0.12_0.1_';
% Input coefficients in DIABLO coordinate

%1) 0.12_10
RI=0.12; Re=1000;Pr=1;
LX=28.28; NX=512;
LY=20;    NY=361;
LZ=7.07;  NZ=128;
%2) 0.12_6
%RI=0.12; Re=1000;Pr=1;
%LX=28.56; NX=512;
%LY=20;    NY=361;
%LZ=7.14;  NZ=128;
%2) 0.12_4
% RI=0.12; Re=1000;Pr=1;
% LX=27.93; NX=512;
% LY=20;    NY=361;
% LZ=6.98;  NZ=128;
%2) 0.12_3
% RI=0.12; Re=1000;Pr=1;
% LX=28.36; NX=512;
% LY=20;    NY=361;
% LZ=7.09;  NZ=128;

% 0.16_3
% RI=0.16; Re=1000;Pr=1;
% LX=28.64; NX=512;
% LY=20;    NY=361;
% LZ=7.16;  NZ=128;
% 0.12_25
% RI=0.12; Re=1000;Pr=1;
% LX=27.93; NX=512;
% LY=20;    NY=361;
% LZ=6.98;  NZ=128;
% 0.12_2
% RI=0.12; Re=1000;Pr=1;
% LX=29.16; NX=512;
% LY=20;    NY=361;
% LZ=7.29;  NZ=128;




%% phase between subharmonic and KH mode at z=0 from 2D slice
x=linspace(0,LX,NX);z=linspace(0,LZ,NZ);
dx=mean(diff(x));
nx=length(x); nz=length(z);

   num=10;
   for i=1:num

filename_mean=[base_dir num2str(i) '/mean.h5'];
filename=[base_dir num2str(i) '/movie.h5'];

file_info=h5info(filename_mean);
att_info=file_info.Groups.Attributes;
nk=att_info.Value;
time=zeros(1,nk);
dphi{i} = zeros(1,nk);
phi_sub{i} = zeros(1,nk);
phi_kh{i} = zeros(1,nk);  

kN = 1/(2*dx); % Nyquist wavenumber

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

 varname=['/time/' timename];            % TIME
 time(k)=h5read(filename_mean,varname);
 end
%if i==2
%k1=1002:1502;
%else i~=2
k1=1:nk;
%end
 for k=1:length(k1)

   if (k1(k)<10)
     timename=['000' int2str(k1(k))];
   elseif (k1(k)<100)
     timename=['00' int2str(k1(k))];
   elseif (k1(k)<1000)
     timename=['0' int2str(k1(k))];
   else
     timename=[int2str(k1(k))];
   end

 varname=['/v_xz/' timename];
 v_xz=h5read(filename,varname);
 
% v' (vertical velocity perturbation)
vprime = v_xz-mean(v_xz,[1,2]);

% taking the spanwise average of the series then take fft
  vpm = mean(vprime,2);                   % series
        % determine the wavenumber
        if mod(nx,2)==0;                       % if N is even
            kj = 2*pi*(-nx/2:nx/2-1)/(nx*dx);         % frequencies
        else mod(nx,2)==1;                     % if N is odd
            kj = 2*pi(-(nx-1)/2:(nx-1)/2)/(nx*dx);   % frequencies
        end
        % fft
        V = zeros(size(kj));
        V = fftshift(fft(vpm,length(kj)))/nx;   % Fourier coeeficient
        % consider only the right half of the wavenumber
        kj = kj(nx/2+2:nx); 
        V = V(nx/2+2:nx);
% compute phase phi for both 1st(subharmonic) and 2nd(KH) mode
   %1) take the ratio of V_KH and V_sub
        C1 = real(V(1));               % real part
        Q1 = imag(V(1));               % imaginery part
        phi_sub{i}(k) = atan2(Q1,C1)/pi;            % phase in radians
        C2 = real(V(2));               % real part
        Q2 = imag(V(2));               % imaginery part
        phi_kh{i}(k) = atan2(Q2,C2)/pi;            % phase in radians
       
   V211 = V(2)./(V(1).^2);
        C211 = real(V211);               % real part
        Q211 = imag(V211);               % imaginery part
        dphi{i}(k) = atan2(Q211,C211)/pi;            % phase in radians
 end
 Time{i} = time;
   end
  
 save dph12_10_sta dphi phi_sub phi_kh Time



