function [Uimx,Uimy,Nump,DispX,DispY,CW]= run_image_matching_uncertainty(im1,im2,window,res,zpad,Zeromean,X,Y,Uin,Vin)
% Input
%im1, im2: deformed images for processing with deform or original images
%for DWO
%window: window size
%res: Window resolution
%zpad: zeropadding
%zeromean: if images are zeromean
%X,Y: vector grid points
%Uin,Vin, velocity from previous pass for DWO.
%Output
%Uimx: Image matching uncertainty in X direction
%Uimy: Image matching uncertainty in Y direction
%Nump: Image matching number of particles detected
% This code is written by Sayantan Bhattacharya


imClass = 'double';

%convert input parameters
im1=cast(im1,imClass);
im2=cast(im2,imClass);
L=size(im1);

%convert to gridpoint list
X=X(:);
Y=Y(:);
% keyboard;
%preallocate velocity fields and grid format
Nx = window(1);
Ny = window(2);

if nargin <=9
    Uin = zeros(length(X),1,imClass);
    Vin = zeros(length(X),1,imClass);
end

%sets up extended domain size
if zpad~=0
    Sy=2*Ny;
    Sx=2*Nx;
else
    Sy=Ny;
    Sx=Nx;
end

%fftshift indicies
% fftindy = [ceil(Sy/2)+1:Sy 1:ceil(Sy/2)];
% fftindx = [ceil(Sx/2)+1:Sx 1:ceil(Sx/2)];

%window masking filter
sfilt1 = windowmask([Sx Sy],[res(1, 1) res(1, 2)]);
sfilt2 = windowmask([Sx Sy],[res(2, 1) res(2, 2)]);

Uimx=zeros(length(X),1,imClass);
Uimy=zeros(length(X),1,imClass);
Nump=zeros(length(X),1,imClass);
DispX = cell(length(X),1);
DispY = cell(length(X),1);
CW    = cell(length(X),1);


for n=1:length(X)
    
    %apply the second order discrete window offset
    x1 = X(n) - floor(round(Uin(n))/2);
    x2 = X(n) +  ceil(round(Uin(n))/2);
    
    y1 = Y(n) - floor(round(Vin(n))/2);
    y2 = Y(n) +  ceil(round(Vin(n))/2);
    
    xmin1 = x1- ceil(Nx/2)+1;
    xmax1 = x1+floor(Nx/2);
    xmin2 = x2- ceil(Nx/2)+1;
    xmax2 = x2+floor(Nx/2);
    ymin1 = y1- ceil(Ny/2)+1;
    ymax1 = y1+floor(Ny/2);
    ymin2 = y2- ceil(Ny/2)+1;
    ymax2 = y2+floor(Ny/2);
    
    %find the image windows
    zone1 = im1( max([1 ymin1]):min([L(1) ymax1]),max([1 xmin1]):min([L(2) xmax1]));
    zone2 = im2( max([1 ymin2]):min([L(1) ymax2]),max([1 xmin2]):min([L(2) xmax2]));
    
    if size(zone1,1)~=Ny || size(zone1,2)~=Nx
        w1 = zeros(Ny,Nx);
        w1( 1+max([0 1-ymin1]):Ny-max([0 ymax1-L(1)]),1+max([0 1-xmin1]):Nx-max([0 xmax1-L(2)]) ) = zone1;
        zone1 = w1;
    end
    if size(zone2,1)~=Ny || size(zone2,2)~=Nx
        w2 = zeros(Ny,Nx);
        w2( 1+max([0 1-ymin2]):Ny-max([0 ymax2-L(1)]),1+max([0 1-xmin2]):Nx-max([0 xmax2-L(2)]) ) = zone2;
        zone2 = w2;
    end
    
    
    if Zeromean==1
        zone1=zone1-mean(mean(zone1));
        zone2=zone2-mean(mean(zone2));
    end
    
    %apply the image spatial filter
    region1 = (zone1).*sfilt1;
    region2 = (zone2).*sfilt2;
    
    
    [deltax,deltay,Np,dispx,dispy,cw]=image_matching_prana(region1,region2);
    Uimx(n)=deltax;
    Uimy(n)=deltay;
    Nump(n)=Np;
    
    if nargout>3
        DispX{n} = dispx;
        DispY{n} = dispy;
        CW{n}    = cw;
    end
end
end