function [tanalpha,tanbeta]=calculate_stereo_angle(calmat,xg,yg,zg,modeltype)
%calculate_stereo_angle calculates the tangent of the camera angles
%   This function calculates the tangent of the camera angles based on the
%   mapping function gradient.
%
%   [tanalpha,tanbeta] = calculate_stereo_angle(calmat,xg,yg,zg,modeltype)
%   where:
%   xg, yg, zg are the locations of the points in the camera object plane
%   modeltype is the form of the camera calibration function with
%   coefficients equal to calmat = [aXcam aYcam]
% 
%   [tanalpha,tanbeta,dFdx1,dFdx2,dFdx3] = calculate_stereo_angle(...)
%   also returns the camera calibration function gradients



%written by Sayantan Bhattacharya January 2016

% % % Viewing angles for camera 1 calculated at gridpoints xgrid,ygrid;
[rows,cols]=size(xg);

aall=calmat;%[aXcam1 aYcam1];
dFdx1=zeros(rows,cols,2);       % the 3rd dimention corresponds to dFdx1 for (X1,Y1)
dFdx2=zeros(rows,cols,2);
dFdx3=zeros(rows,cols,2);

if modeltype==1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mapping the camera coord. to the World Coord. using 1sr order z
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for gg=1:2
        a=aall(:,gg);
        dFdx1(:,:,gg) = a(2) + 2*a(5)*xg + a(6)*yg + a(8)*zg + 3*a(10)*xg.^2 + ...
            2*a(11)*xg.*yg + a(12)*yg.^2 + 2*a(14)*xg.*zg + a(15)*yg.*zg;
        
        dFdx2(:,:,gg) = a(3) + a(6)*xg + 2*a(7)*yg + a(9)*zg + a(11)*xg.^2 + ...
            2*a(12)*xg.*yg + 3*a(13)*yg.^2 + a(15)*xg.*zg + 2*a(16)*yg.*zg;
        
        dFdx3(:,:,gg) = a(4) + a(8)*xg + a(9)*yg + a(14)*xg.^2 + a(15)*xg.*yg + a(16)*yg.^2;
    end
    
elseif modeltype==2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mapping the camera coord. to the World Coord. using 2nd order z
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for gg=1:2
        a=aall(:,gg);
        dFdx1(:,:,gg) = a(2) + 2*a(5).*xg + a(6)*yg + a(8)*zg + 3*a(11)*xg.^2 + 2*a(12)*xg.*yg + ...
            a(13)*yg.^2 + 2*a(15)*xg.*zg + a(16)*yg.*zg + a(18)*zg.^2;
        
        dFdx2(:,:,gg) = a(3) + a(6)*xg + 2*a(7)*yg + a(9)*zg + a(12)*xg.^2 + 2*a(13)*xg.*yg + ...
            3*a(14)*yg.^2 + a(16)*xg.*zg + 2*a(17)*yg.*zg + a(19)*zg.^2;
        
        dFdx3(:,:,gg) = a(4) + a(8)*xg + a(9)*yg + 2*a(10)*zg + a(15)*xg.^2 + a(16)*xg.*yg + ...
            a(17)*yg.^2 + 2*a(18)*xg.*zg + 2*a(19)*yg.*zg;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mapping the camera coord. to the World Coord. using linear interp between cubic xy planes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for gg=1:2
        a=aall(:,gg);
        dFdx1(:,:,gg) = a(2) + 2*a(5)*xg + a(6)*yg + a(8)*zg + 3*a(10)*xg.^2 + ...
            2*a(11)*xg.*yg + a(12)*yg.^2 + 2*a(14)*xg.*zg + a(15)*yg.*zg + ...
            2*a(17)*xg.*yg.*zg + a(18)*yg.^2.*zg + 3*a(19)*xg.^2.*zg;

        dFdx2(:,:,gg) = a(3) + a(6)*xg + 2*a(7)*yg + a(9)*zg + a(11)*xg.^2 + ...
            2*a(12)*xg.*yg + 3*a(13)*yg.^2 + a(15)*xg.*zg + 2*a(16)*yg.*zg + ...
            a(17)*xg.^2.*zg + 2*a(18)*xg.*yg.*zg + 3*a(20)*yg.^2.*zg;

        dFdx3(:,:,gg) = a(4) + a(8)*xg + a(9)*yg + a(14)*xg.^2 + a(15)*xg.*yg + a(16)*yg.^2 + ...
            a(17)*xg.^2.*yg + a(18)*xg.*yg.^2 + a(19)*xg.^3 + a(20)*yg.^3;
    end
end


tanalpha=(((dFdx3(:,:,2).*dFdx2(:,:,1)) - (dFdx2(:,:,2).*dFdx3(:,:,1)))./((dFdx2(:,:,2).*dFdx1(:,:,1)) - (dFdx1(:,:,2).*dFdx2(:,:,1))));
tanbeta =(((dFdx3(:,:,2).*dFdx1(:,:,1)) - (dFdx1(:,:,2).*dFdx3(:,:,1)))./((dFdx1(:,:,2).*dFdx2(:,:,1)) - (dFdx2(:,:,2).*dFdx1(:,:,1))));
% tanalpha2=(((dFdx3(:,:,4).*dFdx2(:,:,3)) - (dFdx2(:,:,4).*dFdx3(:,:,3)))./((dFdx2(:,:,4).*dFdx1(:,:,3)) - (dFdx1(:,:,4).*dFdx2(:,:,3))));
% tanbeta2=(((dFdx3(:,:,4).*dFdx1(:,:,3)) - (dFdx1(:,:,4).*dFdx3(:,:,3)))./((dFdx1(:,:,4).*dFdx2(:,:,3)) - (dFdx2(:,:,4).*dFdx1(:,:,3))));


end