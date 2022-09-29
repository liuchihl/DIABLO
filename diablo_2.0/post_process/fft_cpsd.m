% 1-sided sample estimate of the cross power spectral density (PSD)
function [Sxx,Syy,Sxy,fj]=fft_cpsd(x,y,dt,win)
% This function uses fft routine to calculate F.T and get the 
% cross power spectral density 

% INPUTS:
% x,y:  time series 
% N:    the observation record length
% dt:   time interval: default is set to 1
% Np:   the data record length after zero padding to the time series

% OUTPUT:
% Sp:    PSD
% fj:   the Fourier frequencies 
% written by Chih-Lun Liu, 

% Data Analysis OC683 (reference to Math Camp's note CHP26) 
% Date: May, 2021 
% ----------------------------------------------------

N = length(y);                  % the observation record length

% default
if ~exist('dt'); dt=1; end
if nargin<4; win='rec'; end

    if mod(N,2)==0;                       % if N is even
        fj = (-N/2:N/2-1)/(N*dt);         % frequencies
    else mod(N,2)==1;                     % if N is odd
        fj = (-(N-1)/2:(N-1)/2)/(N*dt);   % frequencies
    end

%%% set the windows 
% bartlett window (no need to boost for cross spectrum)    
    if any(strcmpi(win,'bartlett'));
       y = y.*bartlett(N);      x = x.*bartlett(N);
% hanning window           
    elseif any(strcmpi(win,'hann'));
       y = y.*hann(N);          x = x.*hann(N);
% rectangle window
    else any(strcmpi(win,'rec'));
       y = y.*rectwin(N);       x = x.*rectwin(N);
    end

    
%%% compute auto and cross PSD(2-sided)
    X = fft(x,length(fj))/N;        % FT of x(n)
    X = fftshift(X);                % shifted X
    Y = fft(y,length(fj))/N;        % FT of y(n)
    Y = fftshift(Y);                % shifted Y
%%% auto PSD
    Sxx = X.*conj(X).*N.*dt;        % power spectral density (2-sided)
    Syy = Y.*conj(Y).*N.*dt;        % power spectral density (2-sided)
%%% cross PSD
    Sxy = conj(X).*Y.*N.*dt;        % power spectral density (2-sided)
    
%%% convert PSD to 1-sided: S_prime = 2*S_hat, excluding fj=0 
    Sxx = 2*Sxx(fj>0);    
    Syy = 2*Syy(fj>0);    
    Sxy = Sxy(fj>0);    % no S prime for the cross spectrum, 
                        % so we have to multiply the coherence square 
                        % by 4 in the script
     fj = fj(fj>0);

end
