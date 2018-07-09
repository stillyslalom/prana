function [MI,Nim1,Nim2,INTS,Diapx,Diapy] = MI_Cal_RPC(G,nAuto1,nAuto2,INTS1,INTS2,Diap1x,Diap1y,Diap2x,Diap2y,Sx,Sy,fftindx,fftindy)

%This function is used to calculate the mutual information between two
%consecutive frames using RPC method
%   region 1 and 2 are the particle image within the interrogation window.
%   Sx and Sy are the correlation window size
%   fftindx and fftindy are fftshift indicies
%   We will first find the maximum intensity from the image, and average
%   particle diameter from particle image autocorrelation. A standard
%   gaussian particle is built based on these two parameter. Then the
%   contritbution of one particle is caluclate by take the auotcorrelatio
%   peak of stardard gaussian particle. The primary peak on correlation
%   plane is a summation of the contribution of all correlated particles.
%   MI is the ratio between contribution of all correlated particles and
%   one particle.

% For RPC only magnitude of the cross-correlation is used.

% Diapx=(1/sqrt(2))*mean([Diap1x Diap2x]);
% Diapy=(1/sqrt(2))*mean([Diap1y Diap2y]);
Diapx=(1/sqrt(2))*(Diap1x*Diap2x)^0.5;
Diapy=(1/sqrt(2))*(Diap1y*Diap2y)^0.5;

%make analytical gaussian particle and calculate
%autocorrelaiton of the standard particle
xco = 1:Sx;xco = repmat(xco,[Sx,1]); % build X axis
yco = 1:Sy;yco = repmat(yco',[1,Sy]); % build Y axis
% INTS = (INTS1+INTS2)/2; % intensity of standard Gaussian particle is the mean (maximum) intensity of two frames
INTS = sqrt(INTS1*INTS2); % intensity of standard Gaussian particle is the mean (maximum) intensity of two frames

% % DiaP = (DiaP1+DiaP2)/2; % diameter of standard Gaussian particle is the mean diameter of two frames
% % fg = INTS*exp(-8*((xco-round(Sx/2)).^2+(yco-round(Sy/2)).^2)/DiaP^2); % build the particle intensity distribution based on 2D Gaussian distribution

fg = INTS*exp(-8*((xco-round(Sx/2)).^2/(Diapx^2)+(yco-round(Sy/2)).^2/(Diapy^2))); % build the particle intensity distribution based on 2D Gaussian distribution

fg = fftn(fg,[Sy,Sx]); %FFT
Pg = fg.*conj(fg); %FFT based correlation
Sp = ifftn(Pg,'symmetric'); % convert to time domain
Sp = Sp(fftindy,fftindx); % Sp is the autocorrelation of the standadrd Gaussian particle

 %use cross correaltion plane to get MI
 Gnorm = G-min(G(:)); %minimum correlation subtraction to elminate background noise effect
 [~,Gind] = max(Gnorm(:)); % find the primary peak location 
 shift_locyg = 1+mod(Gind-1,Sy); % find the Y coordinate of the peak
 shift_locxg = ceil(Gind/Sy); % find the X coordinate of the peak
 GXshift = Sy/2+1-shift_locxg; % find the distance between the peak location to the center of correlation plane in X direction
 GYshift = Sx/2+1-shift_locyg; % find the distance between the peak location to the center of correlation plane in Y direction
 SGnorm = circshift(Gnorm,[GYshift,GXshift]); % shift the peak to the center of the plane
 %keyboard;
MI = max(SGnorm(:))/max(Sp(:)); % MI is the ratio between the contribution of all correlated particle and the contribution of one particle
Nim1=max(nAuto1(:))/max(Sp(:));
Nim2=max(nAuto2(:))/max(Sp(:));



end

