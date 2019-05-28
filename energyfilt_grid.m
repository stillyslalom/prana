function [W]=energyfilt_grid(Nx,Ny,d,a,q)
% --- RPC Spectral Filter Subfunction for MTV grids ---

%     This file is part of prana, an open-source GUI-driven program for
%     calculating velocity fields using PIV or PTV.
%
%     Copyright (C) 2012  Virginia Polytechnic Institute and State
%     University
%
%     Copyright 2019. Triad National Security, LLC. All rights reserved.
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


if numel(d) == 1
    d(2) = d;
end

%convert degrees to radians
th1 = a(1) * pi/180;
th2 = a(2) * pi/180;

%assume no aliasing
if nargin<5
    q = 0;
end

%prototype grid crossing
x = -floor(Nx/2):1:(ceil(Nx/2)-1);
y = -floor(Ny/2):1:(ceil(Ny/2)-1);
[X, Y] = meshgrid(x,y);
x0 = 0;
y0 = 0;
I = ( exp( -8*( ((-(X-x0)*sin(th1)+(Y-y0)*cos(th1))).^2 )./d(1).^2 ) ...
    + exp( -8*( ((-(X-x0)*sin(th2)+(Y-y0)*cos(th2))).^2 )./d(2).^2 ) );

%particle-image spectrum
FI = fft2(I);
Ep = FI.*conj(FI);

%aliased particle-image spectrum
%JJC: We haven't derived this, and most of the time q=0 anyway in prana
Ea=zeros(size(Ep));
% % Ea = (pi*255*d^2/8)^2*exp(-d^2*(k1+2*pi).^2/16).*exp(-d^2*(k2+2*pi).^2/16)+...
% %      (pi*255*d^2/8)^2*exp(-d^2*(k1-2*pi).^2/16).*exp(-d^2*(k2+2*pi).^2/16)+...
% %      (pi*255*d^2/8)^2*exp(-d^2*(k1+2*pi).^2/16).*exp(-d^2*(k2-2*pi).^2/16)+...
% %      (pi*255*d^2/8)^2*exp(-d^2*(k1-2*pi).^2/16).*exp(-d^2*(k2-2*pi).^2/16)+...
% %      (pi*255*d^2/8)^2*exp(-d^2*(k1+0*pi).^2/16).*exp(-d^2*(k2+2*pi).^2/16)+...
% %      (pi*255*d^2/8)^2*exp(-d^2*(k1+0*pi).^2/16).*exp(-d^2*(k2-2*pi).^2/16)+...
% %      (pi*255*d^2/8)^2*exp(-d^2*(k1+2*pi).^2/16).*exp(-d^2*(k2+0*pi).^2/16)+...
% %      (pi*255*d^2/8)^2*exp(-d^2*(k1-2*pi).^2/16).*exp(-d^2*(k2+0*pi).^2/16);
% Ea = (pi*255*(d(1)*d(2))/8)^2*exp(-d(2)^2*(k1+2*pi).^2/16).*exp(-d(1)^2*(k2+2*pi).^2/16)+...
%      (pi*255*(d(1)*d(2))/8)^2*exp(-d(2)^2*(k1-2*pi).^2/16).*exp(-d(1)^2*(k2+2*pi).^2/16)+...
%      (pi*255*(d(1)*d(2))/8)^2*exp(-d(2)^2*(k1+2*pi).^2/16).*exp(-d(1)^2*(k2-2*pi).^2/16)+...
%      (pi*255*(d(1)*d(2))/8)^2*exp(-d(2)^2*(k1-2*pi).^2/16).*exp(-d(1)^2*(k2-2*pi).^2/16)+...
%      (pi*255*(d(1)*d(2))/8)^2*exp(-d(2)^2*(k1+0*pi).^2/16).*exp(-d(1)^2*(k2+2*pi).^2/16)+...
%      (pi*255*(d(1)*d(2))/8)^2*exp(-d(2)^2*(k1+0*pi).^2/16).*exp(-d(1)^2*(k2-2*pi).^2/16)+...
%      (pi*255*(d(1)*d(2))/8)^2*exp(-d(2)^2*(k1+2*pi).^2/16).*exp(-d(1)^2*(k2+0*pi).^2/16)+...
%      (pi*255*(d(1)*d(2))/8)^2*exp(-d(2)^2*(k1-2*pi).^2/16).*exp(-d(1)^2*(k2+0*pi).^2/16);

%noise spectrum
En = pi/4*Nx*Ny;

%gridded MTV SNR spectral filter
W  = Ep./((1-q)*En+(q)*Ea);
W  = W/sum(sum(W))*Nx*Ny;


end
