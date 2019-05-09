function [ UX, UY, EPSXX, EPSYY, EPSXY, WXY ] = calcPhaseDisp( p, PHI1_X, PHI1_Y, PHI2_X, PHI2_Y, procedure, maxiter)
% calculate displacement and strain components from phase maps
% 1) input:
% p: grid pitch 
% PHI1_X, PHI1_Y, PHI2_X, PHI2_Y: phase maps in non-deformed (1) and deformed (2) grid images
% procedure: 1 = approximate, 2=PHI2 backdeformed in PHI1 axis, 2=iterative (until relative improvement below "prec" value)
% maxiter: for procedure 2, maximum number of iterations until relative improvement below "stopcrit" value (optional parameter, default: maxiter=1)
% 2) output:
% UX, UY: displacement maps
% EPSXX, EPSYY, EPXY: strain maps
% WXY: rotation map

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

if (nargin==6) 
    maxiter=1; 
end

UX = -( PHI2_X-PHI1_X ) * p/(2*pi);
UY = -( PHI2_Y-PHI1_Y ) * p/(2*pi);

if (procedure==2)
    [Y,X]=ndgrid(1:size(PHI1_X,1),1:size(PHI1_X,2));
    stopcrit=5e-4; % smaller values useless because of interpolation bias
    for n=1:maxiter
        interpPHI2_X= interp2(X,Y,PHI2_X,X+UX,Y+UY,'cubic');
        interpPHI2_Y= interp2(X,Y,PHI2_Y,X+UX,Y+UY,'cubic');
        newUX = -( interpPHI2_X-PHI1_X ) * p/(2*pi);
        newUY = -( interpPHI2_Y-PHI1_Y ) * p/(2*pi);
        deltaX=(newUX(:)-UX(:)); deltaY=(newUY(:)-UY(:));
        UX=newUX;
        UY=newUY;
        if ((norm(deltaX(isfinite(deltaX)))<stopcrit*numel(deltaX))&&(norm(deltaY(isfinite(deltaY)))<stopcrit*numel(deltaY))) 
            break; 
        end
    end
    disp(['procedure 2: ',num2str(n),' iteration(s)  -- delta x:',num2str(norm(deltaX(isfinite(deltaX)))/numel(deltaX)),'  delta y:',num2str(norm(deltaX(isfinite(deltaY)))/numel(deltaY))]); 
end

[GxUX,GyUX]=gradient(UX);
[GxUY,GyUY]=gradient(UY);

EPSXX=GxUX;
EPSYY=GyUY;
EPSXY=(GyUX+GxUY)/2;        
WXY=(GyUX-GxUY)/2; 

end
