function [X,Y,U,V] = PIVgridproc(im1,im2,window,res,zpad,gridSpacing,theta,Zeromean,X,Y,Uin,Vin)
% Process images to find phase

%     This file is part of prana, an open-source GUI-driven program for
%     calculating velocity fields using PIV or PTV.
%
%     Copyright 2018-2019. Triad National Security, LLC. All rights reserved.
%     This program was produced under U.S. Government contract 
%     89233218CNA000001 for Los Alamos National Laboratory (LANL), which is 
%     operated by Triad National Security, LLC for the U.S. Department of 
%     Energy/National Nuclear Security Administration.
%     All rights in the program are reserved by Triad National Security, LLC,
%     and the U.S. Department of Energy/National Nuclear Security Administration. 
%     The Government is granted for itself and others acting on its behalf a 
%     nonexclusive, paid-up, irrevocable worldwide license in this material
%     to reproduce, prepare derivative works, distribute copies to the public, 
%     perform publicly and display publicly, and to permit others to do so.
%     If software is modified to produce derivative works, such modified 
%     software should be clearly marked, so as not to confuse it with the 
%     version available from LANL.
%
%     prana is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

imClass = 'double';

%convert input parameters
%for now, only look at first channel (red) for multicolor images
im1=cast(im1(:,:,1),imClass);
im2=cast(im2(:,:,1),imClass);
L=size(im1);

%Since we process the whole image at once, there really isn't an
%opportunity to use zero mean on a per-window basis.  
%Idea: Can we modify the window function's fft to have a zero component at  
%the zero frequency and get the same effect?
if Zeromean==1
    im1=im1-mean(im1(:));
    im2=im2-mean(im2(:));
end

%convert to gridpoint list
X=X(:);
Y=Y(:);

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

%window masking filter
sfilt1 = buildWindow([res(1, 1) res(1, 2)], [Sx Sy]);
sfilt2 = buildWindow([res(2, 1) res(2, 2)], [Sx Sy]);

% %for consistency, should replace buildWindow with a call to windowmask
% sfilt1 = windowmask([Sx Sy],[res(1, 1) res(1, 2)]);
% sfilt2 = windowmask([Sx Sy],[res(2, 1) res(2, 2)]);
% sfilt1 = sfilt1/sum(sfilt1(:))
% sfilt2 = sfilt2/sum(sfilt2(:))

[pX1, pY1] = phImage(im1, sfilt1, gridSpacing, theta);
[pX2, pY2] = phImage(im2, sfilt2, gridSpacing, theta);

P1_X1 =  pX1;
P1_Y1 =  pY1;
P2_X1 =  pX2;
P2_Y1 =  pY2;

% Unwrap phase
% why cast to single?
P1_X=unwrap_phase(single(P1_X1));
P1_Y=unwrap_phase(single(P1_Y1));
P2_X=unwrap_phase(single(P2_X1));
P2_Y=unwrap_phase(single(P2_Y1));

% Unequvocal correspondence. Needs extension
k=0; l=0;
P2_X=P2_X+l*2*pi;
P2_Y=P2_Y+k*2*pi;

%apply the second order discrete window offset
x1 = X(n) - floor(round(Uin(n))/2);
x2 = X(n) +  ceil(round(Uin(n))/2);

y1 = Y(n) - floor(round(Vin(n))/2);
y2 = Y(n) +  ceil(round(Vin(n))/2);

%clamp x1,x2,y1,y2 to limits of original image
x1(x1<1)    = 1;
y1(y1<1)    = 1;
x1(x1>L(1)) = L(1);
y1(y1>L(2)) = L(2);

x2(x2<1)    = 1;
y2(y2<1)    = 1;
x2(x2>L(1)) = L(1);
y2(y2>L(2)) = L(2);

%truncate phasemap to just the sites we care about
P1_X=P1_X(sub2ind(L,y1,x1));
P1_Y=P1_Y(sub2ind(L,y1,x1));
P2_X=P2_X(sub2ind(L,y2,x2));
P2_Y=P2_Y(sub2ind(L,y2,x2));

U = -( P2_X-P1_X ) * gridSpacing(1)/(2*pi);
V = -( P2_Y-P1_Y ) * gridSpacing(2)/(2*pi);

% [ UX0, UY0, ~, ~, ~, ~ ] = calcPhaseDisp( gridSpacing, P1_X, P1_Y, P2_X, P2_Y, 1, Uin, Vin);

%correct grid angle skew
U=-(U*cosd(theta(1)) + V*cosd(theta(2))); 
V=-(U*sind(theta(1)) + V*sind(theta(2))); 

%add DWO to estimation
U = round(Uin)+U;
V = round(Vin)+V;

end

function [phaseX, phaseY] = phImage(I, g, p, theta)

[X, Y] = meshgrid(0:size(I,2)-1, 0:size(I,1)-1);

% Transform coordinates
X1 = X*cosd(theta(1)) - Y*sind(theta(1));
Y1 = X*cosd(theta(2)) - Y*sind(theta(2)); 
X=X1; Y=Y1;

fc1 = 2*pi/p(1); 
fc2 = 2*pi/p(2); 
% Window FFT
phaseX = winFT(I.* exp(-1i*fc1*X), g, 1);
phaseY = winFT(I.* exp(-1i*fc2*Y), g, 2);

end


function g = buildWindow(R, N)
[XX, YY] = meshgrid(-N(1):N(1),-N(2):N(2));
g = 1/(2*pi*R(1)*R(2))*exp(-((XX.^2)/(2*R(1)^2)+(YY.^2)/(2*R(2)^2)));
g = g/sum(g(:));
end


function phi = winFT(S,g, dirAxis)
% if dirAxis == 2; g=g'; end

% Zero padding
pad = size(S)+size(g)-1;
% FFT along specified direction
Simag = fft2(S, pad(1), pad(2));
% Convolution with Gaussian window
winfft  = Simag .* fft2(g, pad(1), pad(2));
% Inverse FFT
wS = ifft2(winfft);
border=(size(g)-1+mod(size(g)-1,2))/2;
% Remove padding
wS=wS(border(1)+1:border(1)+size(S,1),border(2)+1:border(2)+size(S,2));
% Find phase
phi = angle(wS);
% Find magnitude
mag = abs(wS);
end
