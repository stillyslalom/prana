function[Ixx,Iyy,biasx,biasy,Neff,Autod]=Moment_of_correlation(P21,f1,f2,Sx,Sy,cnorm,D,fftindx,fftindy,G,DXtemp,DYtemp,region1,region2,MIest)
%% This function calculates the standard uncertainty from the cross-correlation plane using Moment of Correlation method.
% written by Sayantan Bhattacharya

% Inputs
%P21= cross-correlation plane in fourier domain
%f1= fft of window1/image1
%f2= fft of window2/image2
%Sx= window size Y
%Sy= window size X
%D=  Approx Diameter of particle image (e.g. [2.8 2.8]), this is used to
%initialize the autocorrelation diameter estimation using subpixel fit

%cnorm= Matrix of ones of size Sx,Sy this is non unity if you want to
%correct for spatial winfow filtering

% fftindx,fftindy= fftshift indices I guess a simple fftshift can be used
% instead, but just to be consistent with prana functions

%G= Cross correlation plane of two windows

%region1, region2= spatially windowed image

%DXtemp,DYtemp= Cross correlation peak diameter estimated from 3pt subpixel
%fit on cross-correlation plane G

%MIest = If MI is already estimated then use that otherwise calculate MI

%Output
%Ixx=PDF diameter in X direction
%Iyy=PDF diameter in X direction
%biasx=bias uncertainty in X direction
%biasy=bias uncertainty in Y direction
%Neff=Effective number of particles
%Autod= Average Autocorrelation diameter


if MIest==-1
    % If Mi has not been estimated
    %Autocorrelations
    P11 = f1.*conj(f1);
    P22 = f2.*conj(f2);
    Auto1 = ifftn(P11,'symmetric');
    Auto2 = ifftn(P22,'symmetric');
    Auto1 = Auto1(fftindy,fftindx);
    Auto2 = Auto2(fftindy,fftindx);
    Auto1=abs(Auto1);
    Auto2=abs(Auto2);
    nAuto1 = Auto1-min(Auto1(:)); % Autocorrelation plane of image 1
    nAuto2 = Auto2-min(Auto2(:)); % Autocorrelation plane of image 2
    
    
    % 3 pt Gaussian fit to Autocorrelation Diameter
    [~,~,~,~,Dauto1x3,Dauto1y3,~]=subpixel(nAuto1,Sx,Sy,cnorm,1,0,D);
    [~,~,~,~,Dauto2x3,Dauto2y3,~]=subpixel(nAuto2,Sx,Sy,cnorm,1,0,D);
    Diap1=sqrt(Dauto1x3*Dauto1y3/2);
    Diap2=sqrt(Dauto2x3*Dauto2y3/2);
    
    %Average Autocorrelation Diameter
    Autod=mean([Diap1 Diap2]);
    
    
    %MI Calculation
    INTS1 = max(region1(:));
    INTS2 = max(region2(:));
    [MI,~,~,~,~,~] = MI_Cal_SCC(G,nAuto1,nAuto2,INTS1,INTS2,Dauto1x3,Dauto1y3,Dauto2x3,Dauto2y3,Sx,Sy,fftindx,fftindy);
else
    % Use estimated MI
    MI=MIest;
    Autod=0;
end
%Cross-correlation subpixel fit
% [U,V,~,~,DXtemp,DYtemp,~]=subpixel(G,Sx,Sy,cnorm,1,0,D);


if isnan(DXtemp)
    DXtemp=sqrt(2)*2.8;
end
if isnan(DYtemp)
    DYtemp=sqrt(2)*2.8;
end

% Estimate the diameter of the cross-correlation peak (with initiation of
% 3pt fit estimates of correlation diameter)
[~,~,~,~,dxc1,dyc1,alpha2]=subpixel(G,Sx,Sy,cnorm,3,0,(1/sqrt(2))*[DXtemp DYtemp]);
% Taking care of peak rotation
DCCx = sqrt( (cos(alpha2)^2*dxc1^2 + sin(alpha2)^2*dyc1^2) );
DCCy = sqrt( (sin(alpha2)^2*dxc1^2 + cos(alpha2)^2*dyc1^2) );
% If estimates are too big than the 3pt fit estimate then default to 3 pt
% fit estimate of the correlation diameter
if DCCx >sqrt(2)*DXtemp
    DCCx=DXtemp;
end
if DCCy >sqrt(2)*DYtemp
    DCCy=DYtemp;
end

%Use Cross-correlation diameter estimate to define convolving Gaussian Diameter
Dconv=(1/sqrt(2))*[DCCx DCCy]; 

% Finding the Phase correlation
W = ones(Sy,Sx);
Wden = sqrt(P21.*conj(P21));
W(Wden~=0) = Wden(Wden~=0);
R = P21./W; % This is effectively the fourier transform of the PDF


% constructing the gaussian which will be convolved with the pdf
[Xt,Yt]=meshgrid(1:Sx,1:Sy);
gfilt=exp(-(4/Dconv(1)^2).*(Xt-Sx/2-1).^2-(4/Dconv(2)^2).*(Yt-Sy/2-1).^2);
Gfilt=abs(fftn(gfilt,[Sy Sx]));

%Convolving the PDF and the Gaussian
G1 = ifftn(R.*Gfilt,'symmetric');
G1 = G1(fftindy,fftindx);
G1 = abs(G1);
G1=G1-min(G1(:));


%subpixel estimation using least squres guassfit for the convolved plane G1
[gpx,gpy,~,~,dx1,dy1,alpha1]=subpixel(G1,Sx,Sy,cnorm,3,0,Dconv);

%find the PDF major and minor axis
Px=real(((dx1(1)^2 - 2*(Dconv(1))^2)).^0.5);
Py=real(((dy1(1)^2 - 2*(Dconv(2))^2)).^0.5);

%If the pdf diameter comes out to be zero or imaginary try 3point fit for
%the convolved plane G1
if Px==0 || Py==0
    % if zero try 3point fit
    [gpx,gpy,~,~,dx1,dy1,~]=subpixel(G1,Sx,Sy,cnorm,1,0,Dconv);
    
    Px=real(((dx1(1)^2 - 2*(Dconv(1))^2)).^0.5);
    Py=real(((dy1(1)^2 - 2*(Dconv(2))^2)).^0.5);
    alpha1=0;
    %project major and minor axis to X and Y axes
    Ixx = sqrt( 1/16 * (cos(alpha1)^2*Px^2 + sin(alpha1)^2*Py^2) );
    Iyy = sqrt( 1/16 * (sin(alpha1)^2*Px^2 + cos(alpha1)^2*Py^2) );
    
    
else
    % If not zero then continue
    %project to axes
    Ixx = sqrt( 1/16 * (cos(alpha1)^2*Px^2 + sin(alpha1)^2*Py^2) );
    Iyy = sqrt( 1/16 * (sin(alpha1)^2*Px^2 + cos(alpha1)^2*Py^2) );
    
end

% The mean particle diameter is estimated from the mean cross-correlation
% peak diameter divided by sqrt(2)
part_dia=(1/sqrt(2))*mean([DCCx DCCy]);

% Effective number of pixels is MI times a circular blob of pixels with
% mean particle diameter
Neff=MI*((pi/4)*(part_dia)^2);

%Bias: This is the peak location of the convolved gaussian plane, typically
%this is done for multipass converged correlation plane in which case the
%peak should be at zero and any deviation is the bias. However this will
%not work for first pass.
biasx=gpx;
biasy=gpy;

% The gradient correction and scaling is done outside PIVwindowed as the
% gradient estimation requires the full velcoity field

% %Gradient correction using velocity gradient field Udiff and Vdiff
% Ixxt= real(sqrt(Ixx.^2 - (Autod.^2/16).*(Udiff).^2 ));
% Iyyt= real(sqrt(Iyy.^2 - (Autod.^2/16).*(Vdiff).^2 ));
% 
% 
% %MC uncertainty after scaling and bias correction
% MCx=sqrt(biasx^2+(Ixxt^2)/Neff);
% MCy=sqrt(biasy^2+(Iyyt^2)/Neff);
end

% JJC: this function is an exact duplicate of the external function.
% Is it necessary?
function [MI,Nim1,Nim2,INTS,Diapx,Diapy] = MI_Cal_SCC(G,nAuto1,nAuto2,INTS1,INTS2,Diap1x,Diap1y,Diap2x,Diap2y,Sx,Sy,fftindx,fftindy)
%[MI,INTS,DiaP] = MI_Cal_SCC(G,nAuto1,nAuto2,INTS1,INTS2,Sx,Sy,fftindx,fftindy)
%[MI,INTS,DiaP] = MI_Cal_SCC( region1,region2,Sx,Sy,fftindx,fftindy)
%This function is used to calculate the mutual information between two
%consecutive frames using SCC method
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

