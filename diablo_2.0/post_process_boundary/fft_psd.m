% 1-sided sample estimate of the power spectral density (PSD)
function [Sp,fj]=fft_psd(y,dt,win,Np)
% This function uses fft routine to calculate F.T and get the 
% power spectral density 

% INPUTS:
% y:    time series 
% dt:   time interval: default is set to 1

% window: different types of filtering
%         default is the rectangle window: 'rec'
%         other options: 
%         'bartlett': Bartlett triangle tapering 
%         'hann': Hanning window

% Np:   the data record length after zero padding to the time series
%       no zero padding in default

% OUTPUT:
% Sp:    1-sided PSD
% fj:    1-sided Fourier frequencies 
% written by Chih-Lun Liu, 

% Data Analysis OC683 (reference to Dudely Chelton's note CHP6) 
% Date: May, 2021 
% ----------------------------------------------------

N = length(y);                  % the observation record length

% default

if ~exist('dt'); dt=1; end
if nargin<3; win='rec'; end
% setting the Fourier frequencies    
    if mod(N,2)==0;                       % if N is even
        fj = (-N/2:N/2-1)/(N*dt);         % frequencies
    else mod(N,2)==1;                     % if N is odd
        fj = (-(N-1)/2:(N-1)/2)/(N*dt);   % frequencies
    end
    
if nargin<4;                              % no zero padding

%%% set the windows 
% bartlett window (boosted)  
    if any(strcmpi(win,'bartlett'));
       var_y = mean(y.^2)-mean(y).^2;   % untapered sample variance
       y = y.*bartlett(N);              % tapered sample with triangle window
       var_yr = mean(y.^2)-mean(y).^2;  % tapered sample variance
       boost_factor = sqrt(var_y./var_yr);

% hanning window           
    elseif any(strcmpi(win,'hann'));
       y=y.*hann(N);

% rectangle window
    else  any(strcmpi(win,'rec'));
       y = y.*rectwin(N);
    end
%%% compute PSD (2-sided)
    Y = fft(y,length(fj))/N;        % FT of y(n)
    PSD = Y.*conj(Y).*N.*dt;        % power spectral density (2-sided)
    PSD = fftshift(PSD);            % shift PSD to [-fN,fN)
%%% convert PSD to 1-sided
    Sp = 2*PSD(find(fj>0):N);    % S_prime = 2*PSD, excluding fj=0 
if any(strcmpi(win,'bartlett'));
    Sp = Sp.*boost_factor.^2;
end
    
    
else             
    % zero padding
    fj = (-Np/2:Np/2-1)./(Np*dt);    % frequencies
    Y = fft(y,length(fj))/Np;       % FT of y(n)
    PSD = Y.*conj(Y).*N.*dt;        % power spectral density
    PSD = fftshift(PSD);              % shift PSD to [-fN,fN)
%%% convert PSD to 1-sided
    Sp = 2*PSD(fj>0);    % S_prime = 2*PSD, excluding fj=0 
end
    fj = fj(fj>0);
end
