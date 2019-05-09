function [U, V] = PIVgridproc(Im1, Im2, gridSpacing, theta, sigG)
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

% Process images to find phase
Im1 = im2double(Im1); Im2 = im2double(Im2);
[pX1, pY1] = phImage(Im1, gridSpacing, theta, sigG);
[pX2, pY2] = phImage(Im2, gridSpacing, theta, sigG);
P1_X1 =  pX1;
P1_Y1 =  pY1;
P2_X1 =  pX2;
P2_Y1 =  pY2;
% Unwrap phase
P1_X=unwrap_phase(single(P1_X1));
P1_Y=unwrap_phase(single(P1_Y1));
P2_X=unwrap_phase(single(P2_X1));
P2_Y=unwrap_phase(single(P2_Y1));

% Unequvocal correspondence. Needs extension
k=0; l=0;
P2_X=P2_X+l*2*pi;
P2_Y=P2_Y+k*2*pi;


[ UX0, UY0, ~, ~, ~, ~ ] = calculate_U_EPS( gridSpacing, P1_X, P1_Y, P2_X, P2_Y, 1);

u=UX0; v=UY0;
U=-(u*cosd(theta(1)) + v*cosd(theta(2))); 
V=-(u*sind(theta(1)) + v*sind(theta(2))); 

end
function [phaseX, phaseY] = phImage(I, p, theta, sigG)

g = buildWindow(p, sigG);

[X, Y] = meshgrid(0:size(I,2)-1, 0:size(I,1)-1);

% Transform coordinates
X1 = X*cosd(theta(1)) - Y*sind(theta(1));
Y1 = X*cosd(theta(2)) - Y*sind(theta(2)); 
X=X1; Y=Y1;

fc = 2*pi/p; 
% Window FFT
phaseX = winFT(I.* exp(-1i*fc*X), g, 1);
phaseY = winFT(I.* exp(-1i*fc*Y), g, 2);

end


function g = buildWindow(gridSpacing, sigG)
[XX, YY] = meshgrid(-sigG:sigG,-sigG:sigG);
g = 1/(2*pi*gridSpacing^2)*exp(-(XX.^2+YY.^2)/(2*gridSpacing^2));
g = g/sum(g(:));
end


function phi = winFT(S,g, dirAxis)
if dirAxis == 2; g=g'; end
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

