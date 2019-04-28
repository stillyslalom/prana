function [Un_alpha1,Un_alpha2,Un_beta1,Un_beta2,tanalpha1,tanalpha2,tanbeta1,tanbeta2]=stereo_angle_uncertainty(caljob,planarjob,dispfield,disp_uncertainty_filedir)
%This function calculates the stereo angle uncertainties 

% Input Variables

%caljob= This is the calibration job structure generated through Prana
%Calibration and conatins the calibration mapping function coefficients and
%the calibration parameters.
%caljob.aXcam1 caljob.aYcam1 caljob.aXcam2 caljob.aYcam2 gives the X and Y
%mapping function coefficents (polynomial model) for camera 1 and 2. This
%can be before or after self-calibration and can be obtained using Prana.

%caljob.modeltype=1 for polynomial model (cubic in x and y and linear in z) 
%caljob.modeltype=2 for polynomial model (cubic in x and y and quadratic in z)

%world coordinates (from prana calibration)
%caljob.allx1data= [x,y,z] world coordinates of calibration grid points for camera1
%caljob.allx2data= [x,y,z] world coordinates of calibration grid points for camera2
%caljob.allX1data= [X,Y] image coordinates of calibration grid points for camera1
%caljob.allX2data= [X,Y] image coordinates of calibration grid points for camera2

%planarjob= The 2D prana jobfile for correlating each camera image pair,
%This argument contains the individual camera imagelist which is used to dewarp the images using the function
%"imagedewarp_alone()". This input can be replaced by "dewarp_grid"
%variable, which essentially contains the dewarped x and y grid on the
%common region of overlap between two cameras.

%dispfield is a structure which contains the disparity field between
%two cameras. The fields in the structure are U,V,X and Y corresponding to
%a disparity vector map obtained using ensemble correlation of two camera
%images at the same time instant.


%disp_uncertainty_filedir= The directory to store the calculated uncertainty in each disparity vector or
%load the uncertainty from a previously calculated matfile in this directory.

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

loadensmdspfld=load([disp_uncertainty_filedir,'ensemble_IM_uncertainty.mat']);
Unwx=loadensmdspfld.Unwx;% Uncertainty in X-direction for disparity vector
Unwy=loadensmdspfld.Unwy;% Uncertainty in Y-direction for disparity vector
dispfield.Unwx=Unwx;
dispfield.Unwy=Unwy;

%% Get Uncertainty in calibration coefficients and the world coordinates xyz 
%% using uncertainty propagation through Triangulation (Both bias and random uncertainty)
%Get dewarped x and y coordinate grid on the overlapping camera domain
[~,dewarp_grid,~]=imagedewarp_alone(caljob,'Willert',planarjob,[1024 1024]);
xgrid1=dewarp_grid.xgrid;
ygrid1=dewarp_grid.ygrid;

%Get z projected locations based on both mean and random uncertainties
%evaluated with ensemble processing
[xg,yg,Uncalcoeff,Ux,Uy,Uz,~,~]=uncertainty_in_xyz_calcoeff(caljob,dispfield,xgrid1,ygrid1);

zg=zeros(size(xg));

%% calculate the uncertainty in the angles using uncertainty in x,y,x and calibration coefficients and the mapping function gradient

% calculate angles
[tanalpha1,tanalpha2,tanbeta1,tanbeta2]=calculate_angle(calmat,xg,yg,zg,caljob.modeltype);

% Find uncertainty in mapping fuction gradient and stereo angles
[Un_alpha1,Un_alpha2,Un_beta1,Un_beta2]=Uncertainty_propagation_through_mapping_function_gradient(calmat,xg,yg,zg,Uncalcoeff,Ux,Uy,Uz);

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
    
end


tanalpha1=(((dFdx3(:,:,2).*dFdx2(:,:,1)) - (dFdx2(:,:,2).*dFdx3(:,:,1)))./((dFdx2(:,:,2).*dFdx1(:,:,1)) - (dFdx1(:,:,2).*dFdx2(:,:,1))));
tanbeta1=(((dFdx3(:,:,2).*dFdx1(:,:,1)) - (dFdx1(:,:,2).*dFdx3(:,:,1)))./((dFdx1(:,:,2).*dFdx2(:,:,1)) - (dFdx2(:,:,2).*dFdx1(:,:,1))));
tanalpha2=(((dFdx3(:,:,4).*dFdx2(:,:,3)) - (dFdx2(:,:,4).*dFdx3(:,:,3)))./((dFdx2(:,:,4).*dFdx1(:,:,3)) - (dFdx1(:,:,4).*dFdx2(:,:,3))));
tanbeta2=(((dFdx3(:,:,4).*dFdx1(:,:,3)) - (dFdx1(:,:,4).*dFdx3(:,:,3)))./((dFdx1(:,:,4).*dFdx2(:,:,3)) - (dFdx2(:,:,4).*dFdx1(:,:,3))));


end

function [Un_alpha1,Un_alpha2,Un_beta1,Un_beta2]=Uncertainty_propagation_through_mapping_function_gradient(calmat,xg,yg,zg,Uncalcoeff,unwx,unwy,unwz)
% This function calculates the uncertainty in the mapping function gradients and subsequesntly the angle uncertainty
%written by Sayantan Bhattacharya January 2016
[r,c]=size(xg);
dFdx1=zeros(r,c,4);       % the 3rd dimention corresponds to dFdx1 for (X1,Y1,X2,Y2)
dFdx2=zeros(r,c,4);
dFdx3=zeros(r,c,4);
Un_dFdx1=zeros(r,c,4);       % the 3rd dimention corresponds to dFdx1 for (X1,Y1,X2,Y2)
Un_dFdx2=zeros(r,c,4);
Un_dFdx3=zeros(r,c,4);

for gg=1:4
    a=calmat(:,gg);
    %        Uc=uncoeff(:,gg);
    
    % Mapping function gradient
    dFdx1(:,:,gg) = a(2) + 2*a(5)*xg + a(6)*yg + a(8)*zg + 3*a(11)*xg.^2 + 2*a(12)*xg.*yg + ...
        a(13)*yg.^2 + 2*a(15)*xg.*zg + a(16)*yg.*zg + a(18)*zg.^2;
    
    dFdx2(:,:,gg) = a(3) + a(6)*xg + 2*a(7)*yg + a(9)*zg + a(12)*xg.^2 + 2*a(13)*xg.*yg + ...
        3*a(14)*yg.^2 + a(16)*xg.*zg + 2*a(17)*yg.*zg + a(19)*zg.^2;
    
    dFdx3(:,:,gg) = a(4) + a(8)*xg + a(9)*yg + 2*a(10)*zg + a(15)*xg.^2 + a(16)*xg.*yg + ...
        a(17)*yg.^2 + 2*a(18)*xg.*zg + 2*a(19)*yg.*zg;
    % Second order derivative
    cfxx=2*a(5) +6*a(11)*xg + 2*a(12)*yg + 2*a(15)*zg ;
    cfxy=a(6) + 2*a(12)*xg +  2*a(13)*yg + a(16)*zg ;
    cfxz=a(8)+ 2*a(15)*xg + a(16)*yg + 2*a(18)*zg;
    
    cfyx=a(6) + 2*a(12)*xg +  2*a(13)*yg + a(16)*zg ;
    cfyy=2*a(7) +2*a(13)*xg + 6*a(14)*yg +  2*a(17)*zg;
    cfyz= a(9)+ a(16)*xg + 2*a(17)*yg + 2*a(19)*zg;
    
    cfzx=a(8)+ 2*a(15)*xg + a(16)*yg + 2*a(18)*zg;
    cfzy=a(9)+ a(16)*xg + 2*a(17)*yg + 2*a(19)*zg;
    cfzz=2*a(10)+ 2*a(18)*xg + 2*a(19)*yg;
    
    %Uncertainty due to Ux, Uy and Uz
    T1dFdx1=sqrt((cfxx.^2).*(unwx.^2)+(cfxy.^2).*(unwy.^2)+(cfxz.^2).*(unwz.^2));
    T1dFdx2=sqrt((cfyx.^2).*(unwx.^2)+(cfyy.^2).*(unwy.^2)+(cfyz.^2).*(unwz.^2));
    T1dFdx3=sqrt((cfzx.^2).*(unwx.^2)+(cfzy.^2).*(unwy.^2)+(cfzz.^2).*(unwz.^2));
    
    
    % Uncertainty due to cal coefficients
    
    a=Uncalcoeff(:,gg); % Uncertainty in mapping function coefficients
        
    T2dFdx1=sqrt(a(2).^2 + (2*a(5)).^2.*xg.^2 + a(6).^2*yg.^2 + a(8).^2*zg.^2 +...
        (3*a(11)).^2*xg.^4 + (2*a(12)).^2*xg.^2.*yg.^2 + a(13).^2*yg.^4 + ...
        (2*a(15)).^2*xg.^2.*zg.^2 + a(16).^2*yg.^2.*zg.^2 + a(18).^2*zg.^4);
    
    T2dFdx2=sqrt(a(3).^2 + a(6).^2*xg.^2 + (2*a(7)).^2*yg.^2 + a(9).^2*zg.^2 +...
        a(12).^2*xg.^4 + (2*a(13)).^2*xg.^2.*yg.^2 + (3*a(14)).^2*yg.^4 + ...
        a(16).^2*xg.^2.*zg.^2 + (2*a(17)).^2*yg.^2.*zg.^2 + a(19).^2*zg.^4);
    
    T2dFdx3=sqrt(a(4).^2 + a(8).^2*xg.^2 + a(9).^2*yg.^2 + (2*a(10)).^2*zg.^2 +...
        a(15).^2*xg.^4 + a(16).^2*xg.^2.*yg.^2 + a(17).^2*yg.^4 + ...
        (2*a(18)).^2*xg.^2.*zg.^2 + (2*a(19)).^2*yg.^2.*zg.^2);
    
    %Total Uncertainty in mapping function gradient
    Un_dFdx1(:,:,gg)=sqrt(T1dFdx1.^2+T2dFdx1.^2);
    Un_dFdx2(:,:,gg)=sqrt(T1dFdx2.^2+T2dFdx2.^2);
    Un_dFdx3(:,:,gg)=sqrt(T1dFdx3.^2+T2dFdx3.^2);

    % Uncertainty in mapping function gradient due to Ux,Uy,Uz
%     Un_dFdx1(:,:,gg)=sqrt(T1dFdx1.^2);
%     Un_dFdx2(:,:,gg)=sqrt(T1dFdx2.^2);
%     Un_dFdx3(:,:,gg)=sqrt(T1dFdx3.^2);
    
    % Uncertainty in mapping function gradient due to Uai's
%     Un_dFdx1(:,:,gg)=sqrt(T2dFdx1.^2);
%     Un_dFdx2(:,:,gg)=sqrt(T2dFdx2.^2);
%     Un_dFdx3(:,:,gg)=sqrt(T2dFdx3.^2);
    
end
% keyboard;

%ANGLE CALCULATION
N1=((dFdx3(:,:,2).*dFdx2(:,:,1))-(dFdx2(:,:,2).*dFdx3(:,:,1)));
D1=(dFdx2(:,:,2).*dFdx1(:,:,1)-dFdx1(:,:,2).*dFdx2(:,:,1));
alpha1 = N1./D1;
%A1=atan2(N1,D1);

N2=((dFdx3(:,:,4).*dFdx2(:,:,3))-(dFdx2(:,:,4).*dFdx3(:,:,3)));
D2=(dFdx2(:,:,4).*dFdx1(:,:,3)-dFdx1(:,:,4).*dFdx2(:,:,3));
alpha2 = N2./D2;
%A2=atan2(N2,D2);

N3=((dFdx3(:,:,2).*dFdx1(:,:,1))-(dFdx1(:,:,2).*dFdx3(:,:,1)));
D3=(dFdx1(:,:,2).*dFdx2(:,:,1)-dFdx2(:,:,2).*dFdx1(:,:,1));
beta1  = N3./D3;
%B1=atan2(N3,D3);

N4=((dFdx3(:,:,4).*dFdx1(:,:,3))-(dFdx1(:,:,4).*dFdx3(:,:,3)));
D4=(dFdx1(:,:,4).*dFdx2(:,:,3)-dFdx2(:,:,4).*dFdx1(:,:,3));
beta2  = N4./D4;
%B2=atan2(N4,D4);


%UNCERTAINTY IN ANGLES
%uncertainty for alpha1
mul1=(1./(1+alpha1.^2)).^2;
term1=(((N1./(D1.^2)).*dFdx2(:,:,2)).^2).*(Un_dFdx1(:,:,1)).^2;
term2=((dFdx3(:,:,2)./D1 + (N1./D1.^2).*dFdx1(:,:,2)).^2).*(Un_dFdx2(:,:,1)).^2;
term3=((dFdx2(:,:,2)./D1).^2).*(Un_dFdx3(:,:,1)).^2;

term4=(((N1./(D1.^2)).*dFdx2(:,:,1)).^2).*(Un_dFdx1(:,:,2)).^2;
term5=((dFdx3(:,:,1)./D1 + (N1./D1.^2).*dFdx1(:,:,1)).^2).*(Un_dFdx2(:,:,2)).^2;
term6=((dFdx2(:,:,1)./D1).^2).*(Un_dFdx3(:,:,2)).^2;

Un_alpha1=(mul1.*(term1+term2+term3+term4+term5+term6)).^0.5;

%uncertainty for alpha2
mul2=(1./(1+alpha2.^2)).^2;
term11=(((N2./(D2.^2)).*dFdx2(:,:,4)).^2).*(Un_dFdx1(:,:,3)).^2;
term12=((dFdx3(:,:,4)./D2 + (N2./D2.^2).*dFdx1(:,:,4)).^2).*(Un_dFdx2(:,:,3)).^2;
term13=((dFdx2(:,:,4)./D2).^2).*(Un_dFdx3(:,:,3)).^2;

term14=(((N2./(D2.^2)).*dFdx2(:,:,3)).^2).*(Un_dFdx1(:,:,4)).^2;
term15=((dFdx3(:,:,3)./D2 + (N2./D2.^2).*dFdx1(:,:,3)).^2).*(Un_dFdx2(:,:,4)).^2;
term16=((dFdx2(:,:,3)./D2).^2).*(Un_dFdx3(:,:,4)).^2;

Un_alpha2=(mul2.*(term11+term12+term13+term14+term15+term16)).^0.5;

%uncertainty for beta1
mul3=(1./(1+beta1.^2)).^2;
term21=(((N3./(D3.^2)).*dFdx1(:,:,2)).^2).*(Un_dFdx2(:,:,1)).^2;
term22=((dFdx3(:,:,2)./D3 + (N3./D3.^2).*dFdx2(:,:,2)).^2).*(Un_dFdx1(:,:,1)).^2;
term23=((dFdx1(:,:,2)./D3).^2).*(Un_dFdx3(:,:,1)).^2;

term24=(((N3./(D3.^2)).*dFdx1(:,:,1)).^2).*(Un_dFdx2(:,:,2)).^2;
term25=((dFdx3(:,:,1)./D3 + (N3./D3.^2).*dFdx2(:,:,1)).^2).*(Un_dFdx1(:,:,2)).^2;
term26=((dFdx1(:,:,1)./D3).^2).*(Un_dFdx3(:,:,2)).^2;

Un_beta1=(mul3.*(term21+term22+term23+term24+term25+term26)).^0.5;

%uncertainty for beta2
mul4=(1./(1+beta2.^2)).^2;
term131=(((N4./(D4.^2)).*dFdx1(:,:,4)).^2).*(Un_dFdx2(:,:,3)).^2;
term132=((dFdx3(:,:,4)./D4 + (N4./D4.^2).*dFdx2(:,:,4)).^2).*(Un_dFdx1(:,:,3)).^2;
term133=((dFdx1(:,:,4)./D4).^2).*(Un_dFdx3(:,:,3)).^2;

term134=(((N4./(D4.^2)).*dFdx1(:,:,3)).^2).*(Un_dFdx2(:,:,4)).^2;
term135=((dFdx3(:,:,3)./D4 + (N4./D4.^2).*dFdx2(:,:,3)).^2).*(Un_dFdx1(:,:,4)).^2;
term136=((dFdx1(:,:,3)./D4).^2).*(Un_dFdx3(:,:,4)).^2;

Un_beta2=(mul4.*(term131+term132+term133+term134+term135+term136)).^0.5;

end
