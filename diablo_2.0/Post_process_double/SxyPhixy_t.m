% This script is part of the boundary proximity project: it calculates the  phase difference between the KH mode and subharmonic mode at the middle of the shear layer
% see crossspec.m for squared coherence (smoothed)
% In this section, data is being extracted and variables are calculated and
% saved
clear;%clc;

% Define file names and directory

  RI=0.16; Re=1000; Pr=1;
  LX=28.56; NX=576;
  LY=30;    NY=613;
  LZ=7.14;  NZ=144;
  D = 2;
base_dir=['/glade/scratch/liuchihl/temp/diablo/DIABLO-master/diablo_2.0/Large_case/doubleshearlayer/Ri016/R1/3D/D_2_1'];
%% phase between subharmonic and KH mode at z=0 from 2D slice
x=linspace(0,LX,NX);z=linspace(0,LZ,NZ);
dx=mean(diff(x));
nx=length(x); nz=length(z);

   num=1;
   for i=1:num
filename_mean=[base_dir '/mean.h5'];
filename=[base_dir '/movie.h5'];

file_info=h5info(filename_mean);
att_info=file_info.Groups.Attributes;
nk=att_info.Value;
time=zeros(1,nk);
dphi{i} = zeros(1,nk);
%phi_sub{i} = zeros(1,nk);
%phi_kh{i} = zeros(1,nk);  

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
 varname=['/gyf/' timename];             % Y-COORDINATE
 gyf(:)=h5read(filename_mean,varname);
 end


k1 = 1:nk;
while any(diff(time)<0.1)
ind1 = find(diff(time)<0.1);
k1(ind1(1)+1) = [];
time(ind1(1)+1) = [];
end
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

 varname=['/v_xy/' timename];
 v_xy=h5read(filename,varname);
 
% v' (vertical velocity perturbation)
% vprime = v_xy-mean(v_xz,[1,2]);

% taking the spanwise average of the series then take fft
a = 0.6585/2;
	% upper billow
	inx1 = find(abs(gyf-D+a)==min(abs(gyf-D+a)));
	v_up(:,k) = v_xy(:,inx1);
        % lower billow
	inx2 = find(abs(gyf+D-a)==min(abs(gyf+D-a)));
        v_low(:,k) = v_xy(:,inx2);
        % determine the wavenumber
        if mod(nx,2)==0;                       % if N is even
            kj = 2*pi*(-nx/2:nx/2-1)/(nx*dx);         % frequencies
        else mod(nx,2)==1;                     % if N is odd
            kj = 2*pi(-(nx-1)/2:(nx-1)/2)/(nx*dx);   % frequencies
        end
        
	% fft the signals
%	UP = fft(v_up(:,k));
%	LOW = fft(v_low(:,k));
	
	% fundamental frequency detection
%	[~, indup] = max(abs(UP));
%	[~, indlow] = max(abs(LOW));	
	% phase difference estimation
%	PhDiff(k) = angle(LOW(indlow)) - angle(UP(indup));
	Time{i} = time;
	[c(:,k),lags(k,:)] = xcorr(v_low(:,k),v_up(:,k));
	[~,inx] = max(abs(c(:,k)));
	Xdiff(k) = lags(k,inx)*dx;
   end
end 
	for k=1:length(time)
		if Xdiff(k) < 0
		Xdiff(k) = Xdiff(k)+LX/2;
		end 
	end
	save lags_innerwaves.mat Xdiff v_up v_low time
%	save PhaseDiffD2 Time v_up v_low PhDiff


