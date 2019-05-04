function [z1grid]= geometricTriangulation(xgrid,x2grid,ygrid,y2grid,zgrid,Dux,Duy,aXcam1,aYcam1,aXcam2,aYcam2,caldata)

%function for calculating local viewing angles and Triangulation
%input calibration matrices and common grid world coordinates
% added caldat in function input argument
%outputs projected z coordinates.

%called in mainselfcal.m

%drawback does not take into account Duy for triangulation... should be
%modified later on;

%written by Sayantan Bhattacharya on 7/19/2013


[rows,cols]=size(zgrid);
%initializing local viewing angles 
alphatan=zeros(rows,cols,2);betatan=zeros(rows,cols,2);
%while (max(max(abs(dz))))>=th %&& (max(max(dz)))<=0.1

% Now calculate viewing angles for both linear z and quadratic z
% calibration models

% Viewing angles for camera 1 calculated at gridpoints xgrid,ygrid;

 aall=[aXcam1 aYcam1 aXcam2 aYcam2];
    dFdx1=zeros(rows,cols,4);       % the 3rd dimention corresponds to dFdx1 for (X1,Y1,X2,Y2)
    dFdx2=zeros(rows,cols,4);
    dFdx3=zeros(rows,cols,4);

    if caldata.modeltype==1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Mapping the camera coord. to the World Coord. using 1sr order z
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for gg=1:2
            a=aall(:,gg);
            dFdx1(:,:,gg) = a(2) + 2*a(5)*xgrid + a(6)*ygrid + a(8)*zgrid + 3*a(10)*xgrid.^2 + ...
                2*a(11)*xgrid.*ygrid + a(12)*ygrid.^2 + 2*a(14)*xgrid.*zgrid + a(15)*ygrid.*zgrid;
            
            dFdx2(:,:,gg) = a(3) + a(6)*xgrid + 2*a(7)*ygrid + a(9)*zgrid + a(11)*xgrid.^2 + ...
                2*a(12)*xgrid.*ygrid + 3*a(13)*ygrid.^2 + a(15)*xgrid.*zgrid + 2*a(16)*ygrid.*zgrid;
            
            dFdx3(:,:,gg) = a(4) + a(8)*xgrid + a(9)*ygrid + a(14)*xgrid.^2 + a(15)*xgrid.*ygrid + a(16)*ygrid.^2;
        end
        
    elseif caldata.modeltype==2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Mapping the camera coord. to the World Coord. using 2nd order z
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for gg=1:2
            a=aall(:,gg);
            dFdx1(:,:,gg) = a(2) + 2*a(5).*xgrid + a(6)*ygrid + a(8)*zgrid + 3*a(11)*xgrid.^2 + 2*a(12)*xgrid.*ygrid + ...
                a(13)*ygrid.^2 + 2*a(15)*xgrid.*zgrid + a(16)*ygrid.*zgrid + a(18)*zgrid.^2;
            
            dFdx2(:,:,gg) = a(3) + a(6)*xgrid + 2*a(7)*ygrid + a(9)*zgrid + a(12)*xgrid.^2 + 2*a(13)*xgrid.*ygrid + ...
                3*a(14)*ygrid.^2 + a(16)*xgrid.*zgrid + 2*a(17)*ygrid.*zgrid + a(19)*zgrid.^2;
            
            dFdx3(:,:,gg) = a(4) + a(8)*xgrid + a(9)*ygrid + 2*a(10)*zgrid + a(15)*xgrid.^2 + a(16)*xgrid.*ygrid + ...
                a(17)*ygrid.^2 + 2*a(18)*xgrid.*zgrid + 2*a(19)*ygrid.*zgrid;
        end
        
    elseif caldata.modeltype==4
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Mapping the camera coord. to the World Coord. using linear interp between cubic xy planes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for gg=1:2
            a=aall(:,gg);
            dFdx1(:,:,gg) = a(2) + 2*a(5)*xgrid + a(6)*ygrid + a(8)*zgrid + 3*a(10)*xgrid.^2 + ...
                2*a(11)*xgrid.*ygrid + a(12)*ygrid.^2 + 2*a(14)*xgrid.*zgrid + a(15)*ygrid.*zgrid + ...
                2*a(17)*xgrid.*ygrid.*zgrid + a(18)*ygrid.^2.*zgrid + 3*a(19)*xgrid.^2.*zgrid;
            
            dFdx2(:,:,gg) = a(3) + a(6)*xgrid + 2*a(7)*ygrid + a(9)*zgrid + a(11)*xgrid.^2 + ...
                2*a(12)*xgrid.*ygrid + 3*a(13)*ygrid.^2 + a(15)*xgrid.*zgrid + 2*a(16)*ygrid.*zgrid + ...
                a(17)*xgrid.^2.*zgrid + 2*a(18)*xgrid.*ygrid.*zgrid + 3*a(20)*ygrid.^2.*zgrid;
            
            dFdx3(:,:,gg) = a(4) + a(8)*xgrid + a(9)*ygrid + a(14)*xgrid.^2 + a(15)*xgrid.*ygrid + a(16)*ygrid.^2 + ...
                a(17)*xgrid.^2.*ygrid + a(18)*xgrid.*ygrid.^2 + a(19)*xgrid.^3 + a(20)*ygrid.^3;
        end
        
    end

%calculating the viewing angles using formula (7) and (8) from Giordano and Astarita's paper "Spatial resolution of the Stereo PIV technique" (2009)

alphatan(:,:,1)=(((dFdx3(:,:,2).*dFdx2(:,:,1)) - (dFdx2(:,:,2).*dFdx3(:,:,1)))./((dFdx2(:,:,2).*dFdx1(:,:,1)) - (dFdx1(:,:,2).*dFdx2(:,:,1))));
betatan(:,:,1)=(((dFdx3(:,:,2).*dFdx1(:,:,1)) - (dFdx1(:,:,2).*dFdx3(:,:,1)))./((dFdx1(:,:,2).*dFdx2(:,:,1)) - (dFdx2(:,:,2).*dFdx1(:,:,1))));

%aall=[aXcam1 aYcam1 aXcam2 aYcam2];
%    dFdx1=zeros(rows,cols,4);       % the 3rd dimention corresponds to dFdx1 for (X1,Y1,X2,Y2)
%    dFdx2=zeros(rows,cols,4);
%    dFdx3=zeros(rows,cols,4);
%Xx1c1=aXcam1(2) + 2*aXcam1(5).*x1


%Viewing angles for camera 2 calculated at gridpoints x2grid=xgrid+Dux,y2grid=ygrid+Duy
    if caldata.modeltype==1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Mapping the camera coord. to the World Coord. using 1sr order z
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for gg=3:4
            a=aall(:,gg);
            dFdx1(:,:,gg) = a(2) + 2*a(5)*x2grid + a(6)*y2grid + a(8)*zgrid + 3*a(10)*x2grid.^2 + ...
                2*a(11)*x2grid.*y2grid + a(12)*y2grid.^2 + 2*a(14)*x2grid.*zgrid + a(15)*y2grid.*zgrid;
            
            dFdx2(:,:,gg) = a(3) + a(6)*x2grid + 2*a(7)*y2grid + a(9)*zgrid + a(11)*x2grid.^2 + ...
                2*a(12)*x2grid.*y2grid + 3*a(13)*y2grid.^2 + a(15)*x2grid.*zgrid + 2*a(16)*y2grid.*zgrid;
            
            dFdx3(:,:,gg) = a(4) + a(8)*x2grid + a(9)*y2grid + a(14)*x2grid.^2 + a(15)*x2grid.*y2grid + a(16)*y2grid.^2;
        end
        
    elseif caldata.modeltype==2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Mapping the camera coord. to the World Coord. using 2nd order z
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for gg=3:4
            a=aall(:,gg);
            dFdx1(:,:,gg) = a(2) + 2*a(5).*x2grid + a(6)*y2grid + a(8)*zgrid + 3*a(11)*x2grid.^2 + 2*a(12)*x2grid.*y2grid + ...
                a(13)*y2grid.^2 + 2*a(15)*x2grid.*zgrid + a(16)*y2grid.*zgrid + a(18)*zgrid.^2;
            
            dFdx2(:,:,gg) = a(3) + a(6)*x2grid + 2*a(7)*y2grid + a(9)*zgrid + a(12)*x2grid.^2 + 2*a(13)*x2grid.*y2grid + ...
                3*a(14)*y2grid.^2 + a(16)*x2grid.*zgrid + 2*a(17)*y2grid.*zgrid + a(19)*zgrid.^2;
            
            dFdx3(:,:,gg) = a(4) + a(8)*x2grid + a(9)*y2grid + 2*a(10)*zgrid + a(15)*x2grid.^2 + a(16)*x2grid.*y2grid + ...
                a(17)*y2grid.^2 + 2*a(18)*x2grid.*zgrid + 2*a(19)*y2grid.*zgrid;
        end
        
    elseif caldata.modeltype==4
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Mapping the camera coord. to the World Coord. using linear interp between cubic xy planes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for gg=3:4
            a=aall(:,gg);
            dFdx1(:,:,gg) = a(2) + 2*a(5)*x2grid + a(6)*y2grid + a(8)*zgrid + 3*a(10)*x2grid.^2 + ...
                2*a(11)*x2grid.*y2grid + a(12)*y2grid.^2 + 2*a(14)*x2grid.*zgrid + a(15)*y2grid.*zgrid + ...
                2*a(17)*x2grid.*y2grid.*zgrid + a(18)*y2grid.^2.*zgrid + 3*a(19)*x2grid.^2.*zgrid;
            
            dFdx2(:,:,gg) = a(3) + a(6)*x2grid + 2*a(7)*y2grid + a(9)*zgrid + a(11)*x2grid.^2 + ...
                2*a(12)*x2grid.*y2grid + 3*a(13)*y2grid.^2 + a(15)*x2grid.*zgrid + 2*a(16)*y2grid.*zgrid + ...
                a(17)*x2grid.^2.*zgrid + 2*a(18)*x2grid.*y2grid.*zgrid + 3*a(20)*y2grid.^2.*zgrid;
            
            dFdx3(:,:,gg) = a(4) + a(8)*x2grid + a(9)*y2grid + a(14)*x2grid.^2 + a(15)*x2grid.*y2grid + a(16)*y2grid.^2 + ...
                a(17)*x2grid.^2.*y2grid + a(18)*x2grid.*y2grid.^2 + a(19)*x2grid.^3 + a(20)*y2grid.^3;
        end
        
    end

alphatan(:,:,2)=(((dFdx3(:,:,4).*dFdx2(:,:,3)) - (dFdx2(:,:,4).*dFdx3(:,:,3)))./((dFdx2(:,:,4).*dFdx1(:,:,3)) - (dFdx1(:,:,4).*dFdx2(:,:,3))));
betatan(:,:,2)=(((dFdx3(:,:,4).*dFdx1(:,:,3)) - (dFdx1(:,:,4).*dFdx3(:,:,3)))./((dFdx1(:,:,4).*dFdx2(:,:,3)) - (dFdx2(:,:,4).*dFdx1(:,:,3))));
%keyboard;

%Display camera angles for reference
    
%         figure(100); subplot(2,2,1);
%         imagesc(atand(alphatan(:,:,1))); colorbar; %caxis([25 30]);
%         title('Camera 1 Angle \alpha1','FontSize',16)
%         subplot(2,2,2);
%         imagesc(atand(alphatan(:,:,2))); colorbar; %caxis([-30 -25]);
%         title('Camera 2 Angle \alpha2','FontSize',16)
%         subplot(2,2,3);
%         imagesc(atand(betatan(:,:,1))); colorbar; %caxis([-2 2]);(end-(zed+10):end-(zed+5))
%         title('Camera 1 Angle \beta1','FontSize',16)
%         subplot(2,2,4);
%         imagesc(atand(betatan(:,:,2))); colorbar; %caxis([-2 2]);
%         title('Camera 1 Angle \beta1','FontSize',16)
%
    %keyboard;
% Triangulation using formula (1) of that paper

%Triangulate based on direction with largest viewing angle
dzX=(-sign(alphatan(1,1)).*Dux)./(abs(alphatan(:,:,1))+abs(alphatan(:,:,2)));
dzY=(sign(betatan(1,1)).*Duy)./(abs(betatan(:,:,1))+abs(betatan(:,:,2)));

%I think the above breaks when the cameras aren't arranged left-right on the same side
dzX2= Dux./(alphatan(:,:,1)-alphatan(:,:,2));
dzY2=-Duy./( betatan(:,:,1)- betatan(:,:,2));

if max(max(abs(alphatan(:,:,1))-abs(alphatan(:,:,2))))>max(max(abs(betatan(:,:,1))-abs(betatan(:,:,2)))) % Previously a bracket was misplaced
    %dz1=(-sign(alphatan(1,1)).*Dux)./(abs(alphatan(:,:,1))+abs(alphatan(:,:,2)));
    dz1=dzX2;
    fprintf('using Dux\n')
else
    %dz2=Duy./(abs(alphatan(:,:,1))+abs(alphatan(:,:,2)));%(abs(betatan(:,:,1))+abs(betatan(:,:,2)));
    %dz1=(sign(betatan(1,1)).*Duy)./(abs(betatan(:,:,1))+abs(betatan(:,:,2)));
    dz1=dzY2;
    fprintf('using Duy\n')
end

%compute the dz that minimizes the residual disparity for both dx and dy,
%weighted by the camera angles
% The following polynomial in dz the cost function for the minimization
% of dr^2 = dx^2 + dy^2.
% The roots will likely be imaginary, the real part will be the dz_min
% [ (tan(a1) - tan(a2))^2 + (tan(b1) - tan(b2))^2, ...
%    2*dy*(tan(b1) - tan(b2)) - 2*dx*(tan(a1) - tan(a2)), ...
%    dx^2 + dy^2]
A = (alphatan(:,:,1)-alphatan(:,:,2)).^2 + (betatan(:,:,1)-betatan(:,:,2)).^2;
B = 2*Duy.*(betatan(:,:,1)-betatan(:,:,2)) - 2*Dux.*(alphatan(:,:,1)-alphatan(:,:,2));
C = Dux.^2 + Duy.^2;

dz2 = real((-B+sqrt(B.^2-4*A.*C))./(2*A));


% %triangulate based on mean of both reconstructions
% dzX=(-sign(alphatan(1,1)).*Dux)./(abs(alphatan(:,:,1))+abs(alphatan(:,:,2)));
% dzY=(sign(betatan(1,1)).*Duy)./(abs(betatan(:,:,1))+abs(betatan(:,:,2)));
% dz1 = (dzX+dzY)/2;


figure(7)
subplot(2,2,1),imagesc(dz2),colorbar,title('dz2')
% subplot(2,2,1),imagesc(dzX),colorbar,title('dzX')
% subplot(2,2,2),imagesc(dzY),colorbar,title('dzY')
subplot(2,2,3),imagesc(dzX2),colorbar,title('dzX2')
subplot(2,2,4),imagesc(dzY2),colorbar,title('dzY2')

% if (mean(abs(Dux(:))))>(mean(abs(Duy(:))))
%     dz1=(-sign(alphatan(1,1)).*Dux)./(abs(alphatan(:,:,1))+abs(alphatan(:,:,2)));
% elseif (mean(abs(Duy(:))))>(mean(abs(Dux(:))))
%     %dz2=Duy./(abs(alphatan(:,:,1))+abs(alphatan(:,:,2)));%(abs(betatan(:,:,1))+abs(betatan(:,:,2)));
%     dz1=(sign(betatan(1,1)).*Duy)./(abs(betatan(:,:,1))+abs(betatan(:,:,2)));
% else
%     dz1=0;
%     fprintf('\n Disparity Map converged \n');
% end
%dz1
%if max(max(abs(Dux)))>max(max(abs(Duy)))
    
%max(max(abs(dz2)))

zgrid=zgrid+dz2;% the new z grid


% figure(21);surf(atand(alphatan(:,:,1)));
% figure(22);surf(atand(alphatan(:,:,2)));
% figure(23);surf(atand(betatan(:,:,1)));
% figure(24);surf(atand(betatan(:,:,2)));
% % % figure(13);surf(dz1);
% % % figure(14);surf(dz2);
% % % figure(15);surf(zgrid);
%  keyboard;

z1grid=zgrid;
%figure(3);surf(z1grid);title('change in z direction');xlabel('x(mm)');ylabel('y(mm)');
%figure(4);surf(z1grid);

end