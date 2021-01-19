function [Un_alpha1,Un_alpha2,Un_beta1,Un_beta2,tanalpha1,tanalpha2,tanbeta1,tanbeta2]=stereo_angle_uncertainty(caljob,dewarp_grid,dispfield,dispuncfield)
%This function calculates the stereo angle uncertainties 
%
% Input Variables
%
%caljob= This is the calibration job structure generated through Prana
%Calibration and conatins the calibration mapping function coefficients and
%the calibration parameters.
%caljob.aXcam1 caljob.aYcam1 caljob.aXcam2 caljob.aYcam2 gives the X and Y
%mapping function coefficents (polynomial model) for camera 1 and 2. This
%can be before or after self-calibration and can be obtained using Prana.
%
%caljob.modeltype=1 for polynomial model (cubic in x and y and linear in z) 
%caljob.modeltype=2 for polynomial model (cubic in x and y and quadratic in z)
%
%world coordinates (from prana calibration)
%caljob.allx1data= [x,y,z] world coordinates of calibration grid points for camera1
%caljob.allx2data= [x,y,z] world coordinates of calibration grid points for camera2
%caljob.allX1data= [X,Y] image coordinates of calibration grid points for camera1
%caljob.allX2data= [X,Y] image coordinates of calibration grid points for camera2
%
%planarjob= The 2D prana jobfile for correlating each camera image pair,
%This argument contains the individual camera imagelist which is used to dewarp the images using the function
%"imagedewarp_alone()". This input can be replaced by "dewarp_grid"
%variable, which essentially contains the dewarped x and y grid on the
%common region of overlap between two cameras.
%
%dispfield is a structure which contains the disparity field between
%two cameras. The fields in the structure are U,V,X and Y corresponding to
%a disparity vector map obtained using ensemble correlation of two camera
%images at the same time instant.
%
%dispuncfield = A structure that contains the calculated uncertainty in each disparity vector.  
%Uncertainties are stored in fields as Unwx and Unwy with same scaling as dispfield
%
%dewarp_grid = a structure with the dewarped x and y coordinate grid on the overlapping camera domain.
%The fields are xgrid and ygrid.

% Output Variables
% tanalpha1,tanalpha2,tanbeta1,tanbeta2= tangent of the stereo angles alpha
% (angle in x-z plane) and beta (angle in y-z plane) for camera 1 and 2

% Un_alpha1,Un_alpha2,Un_beta1,Un_beta2= the uncertainty in stereo
% angles alpha1, beta1, alpha2, beta2

% written by Sayantan Bhattacharya January 2016

%% Get Existing Calibration Matrix
%This assigns the existing calibration matrix to calmat variable
calmat=[caljob.aXcam1 caljob.aYcam1 caljob.aXcam2 caljob.aYcam2];

%% Uncertainty estimated using ensmeble image matching
% Here you load the estimated uncertainty for each disparity vector.
% This is obtained in an ensemble sense i.e. first the disparity field for
% the existing calibration is obtianed using ensemble correlation of camera
% dewarped image pairs. The disparity field is then used to shift the individual
% camera images onto each other(similar to image matching technique) and
% obtain an error distribution for each particle position. The errors are
% stacked up in a histogram for all particles in a window and over all
% frames for which the disparity vector is evaluated. From the ensemble
% error distribution of particle position mismatch the uncertainty in x and
% y for each disparity vector is determined. This is in essence an ensemble
% image matching technique.

dispfield.Unwx=dispuncfield.Unwx;
dispfield.Unwy=dispuncfield.Unwy;

%% Get Uncertainty in calibration coefficients and the world coordinates xyz 
%% using uncertainty propagation through Triangulation (Both bias and random uncertainty)
%Get dewarped x and y coordinate grid on the overlapping camera domain
%[~,dewarp_grid,~]=imagedewarp_alone(caljob,'Willert',planarjob,[1024 1024]);
xgrid=dewarp_grid.xgrid;
ygrid=dewarp_grid.ygrid;

%Get z projected locations based on both mean and random uncertainties
%evaluated with ensemble processing
[xg,yg,Uncalcoeff,Ux,Uy,Uz,~,~]=stereo_uncertainty_in_xyz_calcoeff(caljob,dispfield,xgrid,ygrid);

zg=zeros(size(xg));

%% calculate the uncertainty in the angles using uncertainty in x,y,x and calibration coefficients and the mapping function gradient

% calculate angles
% [tanalpha1,tanalpha2,tanbeta1,tanbeta2]=calculate_angle(calmat,xg,yg,zg,caljob.modeltype);
[tanalpha1,tanbeta1]=calculate_stereo_angle(calmat(:,1:2),xg,yg,zg,caljob.modeltype);
[tanalpha2,tanbeta2]=calculate_stereo_angle(calmat(:,3:4),xg,yg,zg,caljob.modeltype);

% Find uncertainty in mapping fuction gradient and stereo angles
[Un_alpha1,Un_alpha2,Un_beta1,Un_beta2]=stereo_uncertainty_in_map_func_grad(calmat,xg,yg,zg,Uncalcoeff,Ux,Uy,Uz,caljob.modeltype);

% figure;
% subplot(2,2,1);imagesc(rad2deg(atan(alpha1)));colorbar;title('alpha1');
% subplot(2,2,2);imagesc(rad2deg(atan(alpha2)));colorbar;title('alpha2');
% subplot(2,2,3);imagesc(rad2deg(atan(beta1)));colorbar;title('beta1');
% subplot(2,2,4);imagesc(rad2deg(atan(beta2)));colorbar;title('beta2');
% print(gcf,'-dpng',[outputdir,'angles.png'],'-r300');

% figure;
% subplot(2,2,1);imagesc(rad2deg(Un_alpha1));colorbar;title('Un_{alpha1}');caxis([0 0.8]);
% subplot(2,2,2);imagesc(rad2deg(Un_alpha2));colorbar;title('Un_{alpha2}');caxis([0 0.8]);
% subplot(2,2,3);imagesc(rad2deg(Un_beta1));colorbar;title('Un_{beta1}');caxis([0 0.8]);
% subplot(2,2,4);imagesc(rad2deg(Un_beta2));colorbar;title('Un_{beta2}');caxis([0 0.8]);
% print(gcf,'-dpng',[outputdir,'angle_uncertainty.png'],'-r300');



end

%{
function [tanalpha1,tanalpha2,tanbeta1,tanbeta2]=calculate_angle(calmat,xg,yg,zg,modeltype)
% This function calculates the tangent of the camera angles based on the
% mapping function gradient
%written by Sayantan Bhattacharya January 2016
% % % Viewing angles for camera 1 calculated at gridpoints xgrid,ygrid;
[rows,cols]=size(xg);

aall=calmat;%[aXcam1 aYcam1 aXcam2 aYcam2];
dFdx1=zeros(rows,cols,4);       % the 3rd dimention corresponds to dFdx1 for (X1,Y1,X2,Y2)
dFdx2=zeros(rows,cols,4);
dFdx3=zeros(rows,cols,4);

if modeltype==1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mapping the camera coord. to the World Coord. using 1sr order z
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for gg=1:4
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
    for gg=1:4
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
    for gg=1:4
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


tanalpha1=(((dFdx3(:,:,2).*dFdx2(:,:,1)) - (dFdx2(:,:,2).*dFdx3(:,:,1)))./((dFdx2(:,:,2).*dFdx1(:,:,1)) - (dFdx1(:,:,2).*dFdx2(:,:,1))));
tanbeta1=(((dFdx3(:,:,2).*dFdx1(:,:,1)) - (dFdx1(:,:,2).*dFdx3(:,:,1)))./((dFdx1(:,:,2).*dFdx2(:,:,1)) - (dFdx2(:,:,2).*dFdx1(:,:,1))));
tanalpha2=(((dFdx3(:,:,4).*dFdx2(:,:,3)) - (dFdx2(:,:,4).*dFdx3(:,:,3)))./((dFdx2(:,:,4).*dFdx1(:,:,3)) - (dFdx1(:,:,4).*dFdx2(:,:,3))));
tanbeta2=(((dFdx3(:,:,4).*dFdx1(:,:,3)) - (dFdx1(:,:,4).*dFdx3(:,:,3)))./((dFdx1(:,:,4).*dFdx2(:,:,3)) - (dFdx2(:,:,4).*dFdx1(:,:,3))));


end
%}
