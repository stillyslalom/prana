function [xg,yg,Uncal,Ux,Uy,Uz,z1grid,caldatamod]=uncertainty_in_xyz_calcoeff(caldata,dispfield,xgrid1,ygrid1)
% This function calculates the uncertainty in the world coordinates and the
% uncertainty in the calibration coefficients using uncertainty propagation
% through triangulation
% Inputs
%caldata= calibration job containing calibration information and mapping
%functions
%dispfield= disparity vector field U,V based pixel based coordinate grid
%X,Y
%xgrid1,ygrid1= dewarped common grid in physical units

%Outputs
%xg,yg are coordinate grids in physical units where disparity vector is
%evaluated
%Uncal= Uncertainty in calibration coefficients
%Ux,Uy,Uz= Uncertainty in physical coordiantes x,y,z
%z1grid= mean shift in the triangulated z plane due to any residual
%disparity, a measure of bias uncertainty in z
%caldatamod= modified calibration coefficients taking into account the
%uncertainty in the coefficients

%written by Sayantan Bhattacharya on January 2016

%% Calibration information
%world coordinates
allx1data=caldata.allx1data;
allx2data=caldata.allx2data;
%camera image coordinates
allX1data=caldata.allX1data;
allX2data=caldata.allX2data;

method=caldata.modeltype;
optionsls=caldata.optionsls;

%polynomial fitting coefficients
aXcam1=caldata.aXcam1;
aYcam1=caldata.aYcam1;
aXcam2=caldata.aXcam2;
aYcam2=caldata.aYcam2;

calmat=[aXcam1 aYcam1 aXcam2 aYcam2];

%% DISPARITY FIELD information
Dux=dispfield.U(:,:,1);
Duy=dispfield.V(:,:,1);%X and Y disparity matrices need to change correlation points referenced to pixel corners to be referenced to pixel centers
X3=dispfield.X + 0.5;
Y3=dispfield.Y + 0.5;% correlation X and Y grid points (in pixels)
Unwx=dispfield.Unwx;% Random uncertainty using image matching in pixels
Unwy=dispfield.Unwy;

%calculating the length of the common grid in the object plane on which the
%images are dewarped
xmin=min(min(xgrid1));
xmax=max(max(xgrid1));
ymin=min(min(ygrid1));
ymax=max(max(ygrid1));

[Imax1,Jmax1]=size(xgrid1);
% Scaling factor = object plane grid length in x or y(in physical units)mm/ pixel length of images
scalex=(xmax-xmin)/(Jmax1-1); 
scaley=(ymax-ymin)/(Imax1-1);
% [scalex scaley]
%Do coordinate transform between vectors location reported as image indices 
% (From 1 to Imax1) to world coordinates:
xg = scalex * (X3-1) + xmin; %(in physical units or mm)
yg = scaley * (Y3-1) + ymin;
zgrid=zeros(size(xg));
%% For Bias uncertainty z location
%scaling the disparity to physical dimensions 
Duxb=Dux.*(scalex);
Duyb=Duy.*(scaley);

% Removing outliers at the edges
stdx=rms(Duxb(:)-mean(Duxb(:)));
stdy=rms(Duyb(:)-mean(Duyb(:)));
meanx=mean(Duxb(:));
meany=mean(Duyb(:));

for i=1:size(xg,1)
    for j=1:size(xg,2)
        
        if abs(Duxb(i,j))<abs(meanx)-2*stdx
            Duxb(i,j)=meanx;
        end
        if abs(Duyb(i,j))<abs(meany)-2*stdy
            Duyb(i,j)=meany;
        end
    end
end


%world grid points shifted by the amount of disparity to get the
%locations at which camera 2 local viewing angle are calculated.
x1grid=xg-Duxb./2;
x2grid=xg+Duxb./2;
y1grid=yg-Duyb./2;
y2grid=yg+Duyb./2; 

%Doing Triangulation of bias uncertainty term
[z1grid]= geometricTriangulation(x1grid,x2grid,y1grid,y2grid,zgrid,Duxb,Duyb,aXcam1,aYcam1,aXcam2,aYcam2,caldata); %ouput is the projected z world points

zbias=z1grid(:); % z bias  uncertainty
%% For Random uncertainty z location
%scaling the disparity to physical dimensions 
Duxr=Unwx.*(scalex);
Duyr=Unwy.*(scaley);
%world grid points shifted by the amount of disparity to get the
%locations at which camera 2 local viewing angle are calculated.
x3grid=xg-Duxr./2;
x4grid=xg+Duxr./2;
y3grid=yg-Duyr./2;
y4grid=yg+Duyr./2; 

%Doing Triangulation of random uncertainty term
[Unwz]= geometricTriangulation(x3grid,x4grid,y3grid,y4grid,zgrid,Duxr,Duyr,aXcam1,aYcam1,aXcam2,aYcam2,caldata); %ouput is the projected z world points

zrand=abs(Unwz(:)); % z rand  uncertainty
%% uncertainty in x, y and z world coordinates

Ux=sqrt(Duxb.^2+Duxr.^2);
Uy=sqrt(Duyb.^2+Duyr.^2);
Uz1=sqrt(zbias.^2+zrand.^2); % Z coordinate uncertainty
Uz=reshape(Uz1,size(xg)); %Reshaping the Z coordinate uncertainty to a matrix form

%% new z plane for bias uncertainty
% If there is a mean shift in the z-plane because of the residual
% uncertainty then evaluate the corresponding bias uncertainty in
% calibration coefficients

xg1=xg(:);
yg1=yg(:);
zg1=zbias;

a1=length(xg1);a2=sum(xg1);a3=sum(yg1);
a4=sum(xg1.*yg1); a7=sum(xg1.^2);a8=sum(yg1.^2);
a5=sum(xg1.*zg1);a6=sum(yg1.*zg1); a9=sum(zg1);

AA=[a7 a4 a2; a4 a8 a3; a2 a3 a1];
bb=[a5;a6;a9];

% New z plane coefficients
coeff=AA\bb; 

% New X-Z plane
% newzx=coeff(1).*xg1 + coeff(3);

Ndof=length(xg1)-length(coeff);

% Error in fit
errfit=zbias-(coeff(1).*xg1+coeff(2).*yg1+coeff(3));
% RMS error in fit or fit uncertainty
zfit=sqrt(sum(errfit.^2)/Ndof);
zstd1=zfit.*ones(size(xg1));

% Total random uncertainty in z as combination of random uncertainty due to
% uncertainty in disparity vector and the fit uncertainty
zrand2=sqrt(zrand(:).^2+zstd1.^2);


%% Plane fitting coefficients and its uncertainty
C=inv(AA);

Uncoeff1=sqrt(sum(  (C(1,1)*xg1+C(1,2)*yg1+C(1,3)).^2.*(zrand2.^2)    ));
Uncoeff2=sqrt(sum(  (C(2,1)*xg1+C(2,2)*yg1+C(2,3)).^2.*(zrand2.^2)    ));
Uncoeff3=sqrt(sum(  (C(3,1)*xg1+C(3,2)*yg1+C(3,3)).^2.*(zrand2.^2)    ));

% positive and negative z planes about the mean 
coeffp=coeff-[Uncoeff1;Uncoeff2;Uncoeff3];
coeffn=coeff+[Uncoeff1;Uncoeff2;Uncoeff3];

%%%%%%%%%
%{
% % zp=zeros(size(zrand));
% % zn=zeros(size(zrand));
% % for lv=1:length(xg1)
% %     
% %     if abs(zg1(lv))>=abs(newzx(lv))
% %         zp(lv)=zbias(lv)+zrand(lv);
% %         zn(lv)=zbias(lv)-zrand(lv);
% %     elseif abs(zg1(lv))<abs(newzx(lv))
% %         zp(lv)=zbias(lv)-zrand(lv);
% %         zn(lv)=zbias(lv)+zrand(lv);
% %     end
% % end


% zp=zbias+zrand;
% zn=zbias-zrand;

%fitting a plane to the projected z points for random uncertainty

% zg1p=zp;
% zg1n=zn;

% % a5p=sum(xg1.*zg1p);a6p=sum(yg1.*zg1p); a9p=sum(zg1p);
% % a5n=sum(xg1.*zg1n);a6n=sum(yg1.*zg1n); a9n=sum(zg1n);
% % 
% % bbp=[a5p;a6p;a9p];
% % bbn=[a5n;a6n;a9n];
% % 
% % % coefficients for the plane Z= Ax+By+C
% % coeffp=AA\bbp;
% % coeffn=AA\bbn;

%%%%%%%%%%%%%
% % newzxp=coeffp(1).*xg1 + coeffp(3);
% % newzxn=coeffn(1).*xg1 + coeffn(3);

%% plot z locations

% % ms1=8;
% % ms2=10;
% % skp=160;
% % skp2=480;
% % figure;hold on;
% % plot(xg1(1:skp:end),zbias(1:skp:end),'k.','markersize',ms1);
% % plot(xg1(1:skp:end),zp(1:skp:end),'r.','markersize',ms1);
% % plot(xg1(1:skp:end),zn(1:skp:end),'r.','markersize',ms1);
% % plot(xg1(1:skp2:end),newzx(1:skp2:end),'k*','markersize',5);
% % plot(xg1(1:skp2:end),newzxp(1:skp2:end),'b.','markersize',ms2);
% % plot(xg1(1:skp2:end),newzxn(1:skp2:end),'b.','markersize',ms2);
% % hold off;
% % % axis([-inf inf -3.2 -2.9]);
% % title({'World coordinate system least squares fit','and its uncertainty'});
% % hltag=legend('z_b','z_b+z_r','z_b-z_r','lsqfit z^{wc}_{b}','lsqfit z^{wc}_{b}+z^{wc}_{r}','lsqfit z^{wc}_{b}-z^{wc}_{r}');
% % set(hltag,'location','eastoutside');
% % set(hltag,'FontSize',12);
% % xlabel('x grid points');ylabel('z position (mm)');

% axis([-inf inf -3.07 -3.00]);
% keyboard;
% print(gcf,'-dpng','calzuncertainty_01_24_withselfcal.png','-r300');
% print(gcf,'-dpng','calzuncertainty_01_24_noselfcal.png','-r300');
% print(gcf,'-dpng','calzuncertainty_49_72_withoffset.png','-r300');
%}
%%%%%%%%%
%% Roation and translation of world coordinates of calibration grid points based on New Z plane locations (mean location +- location due to random uncertainty in z)
alpha=atan(coeff(1));beta=atan(coeff(2));                                  % x-z and y-z plane slopes
Roty=[cos(alpha) 0 -sin(alpha);0 1 0;sin(alpha) 0 cos(alpha)];             %rotation about x and y axis
Rotx=[1 0 0;0 cos(beta) -sin(beta);0 sin(beta) cos(beta)];

alphap=atan(coeffp(1));betap=atan(coeffp(2));                                  % x-z and y-z plane slopes
Rotyp=[cos(alphap) 0 -sin(alphap);0 1 0;sin(alphap) 0 cos(alphap)];             %rotation about x and y axis
Rotxp=[1 0 0;0 cos(betap) -sin(betap);0 sin(betap) cos(betap)];

alphan=atan(coeffn(1));betan=atan(coeffn(2));                                  % x-z and y-z plane slopes
Rotyn=[cos(alphan) 0 -sin(alphan);0 1 0;sin(alphan) 0 cos(alphan)];             %rotation about x and y axis
Rotxn=[1 0 0;0 cos(betan) -sin(betan);0 sin(betan) cos(betan)];

%translation along z
tz=[0;0;coeff(3)]; 
tzp=[0;0;coeffp(3)]; 
tzn=[0;0;coeffn(3)]; 

l1=length(allx1data(:,1));
l2=length(allx2data(:,1));

%Modifying the calibration plane world coordinates  for camera1 and camera2
%such that the new fitted z plane becomes z=0 or in other words the
%original calibration plane world z coordinates are calculated with respect
%to the new fitted plane.
ztrans1=Roty'*Rotx'*[allx1data(:,1)';allx1data(:,2)';allx1data(:,3)'] - [tz(1).*ones(1,l1);tz(2).*ones(1,l1);tz(3).*ones(1,l1)];
ztrans2=Roty'*Rotx'*[allx2data(:,1)';allx2data(:,2)';allx2data(:,3)'] - [tz(1).*ones(1,l2);tz(2).*ones(1,l2);tz(3).*ones(1,l2)];

ztrans1p=Rotyp'*Rotxp'*[allx1data(:,1)';allx1data(:,2)';allx1data(:,3)'] - [tzp(1).*ones(1,l1);tz(2).*ones(1,l1);tz(3).*ones(1,l1)];
ztrans2p=Rotyp'*Rotxp'*[allx2data(:,1)';allx2data(:,2)';allx2data(:,3)'] - [tzp(1).*ones(1,l2);tz(2).*ones(1,l2);tz(3).*ones(1,l2)];

ztrans1n=Rotyn'*Rotxn'*[allx1data(:,1)';allx1data(:,2)';allx1data(:,3)'] - [tzn(1).*ones(1,l1);tz(2).*ones(1,l1);tz(3).*ones(1,l1)];
ztrans2n=Rotyn'*Rotxn'*[allx2data(:,1)';allx2data(:,2)';allx2data(:,3)'] - [tzn(1).*ones(1,l2);tz(2).*ones(1,l2);tz(3).*ones(1,l2)];


%% Find cal coefficients for new z planes
caldatamod.allx1data=ztrans1';
caldatamod.allx2data=ztrans2'; %outputs the modified planes
caldatamod.allx1datap=ztrans1p';
caldatamod.allx2datap=ztrans2p';
caldatamod.allx1datan=ztrans1n';
caldatamod.allx2datan=ztrans2n';
%calculate new polynomial transform coefficients
[~,~,TaXcam1, TaYcam1, TaXcam2, TaYcam2,~,Uncalfit1]=fitmodels_mod(caldatamod.allx1data,...
    caldatamod.allx2data,allX1data,allX2data,method,optionsls);
caldatamod.aXcam1=TaXcam1;
caldatamod.aYcam1=TaYcam1;
caldatamod.aXcam2=TaXcam2;
caldatamod.aYcam2=TaYcam2;

% keyboard;
% % [calmat2]=calibration_lsfit(caldatamod.allx1data,caldatamod.allx2data,allX1data,allX2data,method);
% [calmat1]=calibration_lsfit(caldatamod.allx1data,allX1data,method);
% [calmat2]=calibration_lsfit(caldatamod.allx2data,allX2data,method);



[~,~,TaXcam1, TaYcam1, TaXcam2, TaYcam2,~,~]=fitmodels_mod(caldatamod.allx1datap,...
    caldatamod.allx2datap,allX1data,allX2data,method,optionsls);
caldatamod.aXcam1p=TaXcam1;
caldatamod.aYcam1p=TaYcam1;
caldatamod.aXcam2p=TaXcam2;
caldatamod.aYcam2p=TaYcam2;

[~,~,TaXcam1, TaYcam1, TaXcam2, TaYcam2,~,~]=fitmodels_mod(caldatamod.allx1datan,...
    caldatamod.allx2datan,allX1data,allX2data,method,optionsls);
caldatamod.aXcam1n=TaXcam1;
caldatamod.aYcam1n=TaYcam1;
caldatamod.aXcam2n=TaXcam2;
caldatamod.aYcam2n=TaYcam2;

%mapping function coefficients due to mean shift in z plane location
biascal=[caldatamod.aXcam1 caldatamod.aYcam1 caldatamod.aXcam2 caldatamod.aYcam2];
%mapping function coefficients due to mean + positive sigma shift in z plane location
randpcal=[caldatamod.aXcam1p caldatamod.aYcam1p caldatamod.aXcam2p caldatamod.aYcam2p];
%mapping function coefficients due to mean - positive sigma shift in z plane location
randncal=[caldatamod.aXcam1n caldatamod.aYcam1n caldatamod.aXcam2n caldatamod.aYcam2n];

%Bias uncertianty
Unbiascal=abs(calmat-biascal);
%Random uncertainty
Unrandcal=abs(randpcal-randncal)./2;
%Total Uncertainty in calibration mapping function coefficients
Uncal=sqrt(Unbiascal.^2+Unrandcal.^2+Uncalfit1.^2);
end

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
        for gg=1:4
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
        for gg=1:4
            a=aall(:,gg);
            dFdx1(:,:,gg) = a(2) + 2*a(5).*xgrid + a(6)*ygrid + a(8)*zgrid + 3*a(11)*xgrid.^2 + 2*a(12)*xgrid.*ygrid + ...
                a(13)*ygrid.^2 + 2*a(15)*xgrid.*zgrid + a(16)*ygrid.*zgrid + a(18)*zgrid.^2;
            
            dFdx2(:,:,gg) = a(3) + a(6)*xgrid + 2*a(7)*ygrid + a(9)*zgrid + a(12)*xgrid.^2 + 2*a(13)*xgrid.*ygrid + ...
                3*a(14)*ygrid.^2 + a(16)*xgrid.*zgrid + 2*a(17)*ygrid.*zgrid + a(19)*zgrid.^2;
            
            dFdx3(:,:,gg) = a(4) + a(8)*xgrid + a(9)*ygrid + 2*a(10)*zgrid + a(15)*xgrid.^2 + a(16)*xgrid.*ygrid + ...
                a(17)*ygrid.^2 + 2*a(18)*xgrid.*zgrid + 2*a(19)*ygrid.*zgrid;
        end
        
    end

%calculating the viewing angles using formula (7) and (8) from Giordano and Astarita's paper "Spatial resolution of the Stereo PIV technique" (2009)

alphatan(:,:,1)=(((dFdx3(:,:,2).*dFdx2(:,:,1)) - (dFdx2(:,:,2).*dFdx3(:,:,1)))./((dFdx2(:,:,2).*dFdx1(:,:,1)) - (dFdx1(:,:,2).*dFdx2(:,:,1))));
betatan(:,:,1)=(((dFdx3(:,:,2).*dFdx1(:,:,1)) - (dFdx1(:,:,2).*dFdx3(:,:,1)))./((dFdx1(:,:,2).*dFdx2(:,:,1)) - (dFdx2(:,:,2).*dFdx1(:,:,1))));

%aall=[aXcam1 aYcam1 aXcam2 aYcam2];
    dFdx1=zeros(rows,cols,4);       % the 3rd dimention corresponds to dFdx1 for (X1,Y1,X2,Y2)
    dFdx2=zeros(rows,cols,4);
    dFdx3=zeros(rows,cols,4);
%Xx1c1=aXcam1(2) + 2*aXcam1(5).*x1


%Viewing angles for camera 2 calculated at gridpoints x2grid=xgrid+Dux,y2grid=ygrid+Duy
    if caldata.modeltype==1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Mapping the camera coord. to the World Coord. using 1sr order z
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for gg=1:4
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
        for gg=1:4
            a=aall(:,gg);
            dFdx1(:,:,gg) = a(2) + 2*a(5).*x2grid + a(6)*y2grid + a(8)*zgrid + 3*a(11)*x2grid.^2 + 2*a(12)*x2grid.*y2grid + ...
                a(13)*y2grid.^2 + 2*a(15)*x2grid.*zgrid + a(16)*y2grid.*zgrid + a(18)*zgrid.^2;
            
            dFdx2(:,:,gg) = a(3) + a(6)*x2grid + 2*a(7)*y2grid + a(9)*zgrid + a(12)*x2grid.^2 + 2*a(13)*x2grid.*y2grid + ...
                3*a(14)*y2grid.^2 + a(16)*x2grid.*zgrid + 2*a(17)*y2grid.*zgrid + a(19)*zgrid.^2;
            
            dFdx3(:,:,gg) = a(4) + a(8)*x2grid + a(9)*y2grid + 2*a(10)*zgrid + a(15)*x2grid.^2 + a(16)*x2grid.*y2grid + ...
                a(17)*y2grid.^2 + 2*a(18)*x2grid.*zgrid + 2*a(19)*y2grid.*zgrid;
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
%     keyboard;

alpha1=atand(alphatan(:,:,1));
alpha2=atand(alphatan(:,:,2));
beta1=atand(betatan(:,:,1));
beta2=atand(betatan(:,:,2));

% Triangulation using formula (1) of that paper

if max(abs(alpha1(:))+abs(alpha2(:))) > max(abs(beta1(:))+abs(beta2(:))) 
    % if cameras are placed along x axis when alpha>beta
    if mean(alpha1(:))<mean(alpha2(:)) 
        %1
        dz1=-(Dux)./(abs(alphatan(:,:,1))+abs(alphatan(:,:,2)));
    elseif mean(alpha1(:))>mean(alpha2(:)) 
        %2
        dz1=(Dux)./(abs(alphatan(:,:,1))+abs(alphatan(:,:,2)));
    end
    
elseif max(abs(alpha1(:))+abs(alpha2(:))) < max(abs(beta1(:))+abs(beta2(:)))   
    % if cameras are placed along y axis when alpha<beta
    if mean(beta1(:))<mean(beta2(:)) 
        %3
        dz1=-(Duy)./(abs(betatan(:,:,1))+abs(betatan(:,:,2)));
    elseif mean(beta1(:))>mean(beta2(:)) 
        %4
        dz1=(Duy)./(abs(betatan(:,:,1))+abs(betatan(:,:,2)));
    end
    
end


zgrid=zgrid+dz1;% the new z grid

z1grid=zgrid;
%figure(3);surf(z1grid);title('change in z direction');xlabel('x(mm)');ylabel('y(mm)');
%figure(4);surf(z1grid);

end

function [a_cam1, a_cam2, aXcam1, aYcam1, aXcam2, aYcam2, convergemessage,Uncal] = fitmodels_mod(allx1data,...
    allx2data,allX1data,allX2data,modeltype,~)
% function [a_cam1 a_cam2 aXcam1 aYcam1 aXcam2 aYcam2 convergemessage]=fitcameramodels(allx1data,...
%     allx2data,allX1data,allX2data,modeltype,optionsls)
% This function ....
% Inputs:
%   allx1data:          Location of target markers (x,y,z) in world coord.
%                       for camera 1.
%   allx2data:          Location of target markers (x,y,z) in world coord.
%                       for camera 2.
%   allX1data:          Location of target markers (x,y) in pixel coord.
%                       for camera 1.
%   allX2data:          Location of target markers (x,y) in pixel coord.
%                       for camera 2.
%   Modeltype:          A switch to choose which camera model that will be
%                         used for the 3D fit.  1 = Cubic XY, Linear Z, 2 =
%                         Cubic XY, Quadratic Z, & 3 = Camera Pinhole Model
%                         (DLT).
%   optionsls:          These are the options for the least squares solver.
%                         This variable is determined by optimset
% 
%  Outputs:
%    a_cam1:            Polinomial coeff for camera 1 using the pinhole 
%                         camera model.
%    a_cam2:            Polinomial coeff for camera 2 using the pinhole 
%                         camera model.
%    a_Xcam1:           Polinomial coeff for the X component of camera 1
%                         using the cubic mapping function.
%    a_Ycam1:           Polinomial coeff for the Y component of camera 1
%                         using the cubic mapping function.
%    a_Xcam2:           Polinomial coeff for the X component of camera 2
%                         using the cubic mapping function.
%    a_Ycam2:           Polinomial coeff for the Y component of camera 2
%                         using the cubic mapping function.
%    convergemessage:   Provides a message about the convergence of the 
%                         non-linear least squares solver.  This variable 
%                         has 4 parts, (X and Y for both cameras).  The tag
%                         following the Yes or No is from the lsqnonlin
%                         function built into matlab.  The next to parts
%                         are RMS in pixels and the R value.

% Writen by M. Brady
% Edited and Commented by S. Raben

%     This file is part of prana, an open-source GUI-driven program for
%     calculating velocity fields using PIV or PTV.
%
%     Copyright (C) 2014  Virginia Polytechnic Institute and State
%     University
% 
%     Copyright 2014.  Los Alamos National Security, LLC. This material was
%     produced under U.S. Government contract DE-AC52-06NA25396 for Los 
%     Alamos National Laboratory (LANL), which is operated by Los Alamos 
%     National Security, LLC for the U.S. Department of Energy. The U.S. 
%     Government has rights to use, reproduce, and distribute this software.
%     NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY
%     WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF
%     THIS SOFTWARE.  If software is modified to produce derivative works,
%     such modified software should be clearly marked, so as not to confuse
%     it with the version available from LANL.
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


% save('allx1dataFLS.mat','allx1data');
% save('allx2dataFLS.mat','allx2data');
% save('allX1dataFLS.mat','allX1data');
% save('allX2dataFLS.mat','allX2data');

% save('allx1datanew.mat','allx1data');
% save('allx2datanew.mat','allx2data');
% save('allX1datanew.mat','allX1data');
% save('allX2datanew.mat','allX2data');

% save('allx1datagood.mat','allx1data');
% save('allx2datagood.mat','allx2data');
% save('allX1datagood.mat','allX1data');
% save('allX2datagood.mat','allX2data');

%keyboard;
%  allx1data(:,1)=(allx1data(:,1)-mean(allx1data(:,1)))./(1*max(abs(allx1data(:,1))));
%  allx1data(:,2)=(allx1data(:,2)-mean(allx1data(:,2)))./(1*max(abs(allx1data(:,2))));
%  allx1data(:,3)=(allx1data(:,3)-mean(allx1data(:,3)))./(1*max(abs(allx1data(:,3))));
%  
%  allx2data(:,1)=(allx2data(:,1)-mean(allx2data(:,1)))./(1*max(abs(allx2data(:,1))));
%  allx2data(:,2)=(allx2data(:,2)-mean(allx2data(:,2)))./(1*max(abs(allx2data(:,2))));
%  allx2data(:,3)=(allx2data(:,3)-mean(allx2data(:,3)))./(1*max(abs(allx2data(:,3))));
%  
%  allX1data(:,1)=(allX1data(:,1)-mean(allX1data(:,1)))./(1*max(abs(allX1data(:,1))));
%  allX1data(:,2)=(allX1data(:,2)-mean(allX1data(:,2)))./(1*max(abs(allX1data(:,2))));
%   
%  allX2data(:,1)=(allX2data(:,1)-mean(allX2data(:,1)))./(1*max(abs(allX2data(:,1))));
%  allX2data(:,2)=(allX2data(:,2)-mean(allX2data(:,2)))./(1*max(abs(allX2data(:,2))));

 
 
 
optionsls=optimset('MaxIter',30000,'MaxFunEvals',30000,'TolX',1e-11,'TolFun',1e-7,...
        'LargeScale','off','Display','off','Algorithm','levenberg-marquardt');


% Matlab perfers 'Levenberg-Marquardt' over 'levenbergmarquardt' but the
% optimset function hasn't been updated yet.  This should be fixed in later
% versions of matlab
% fprintf('LevenbergMarquart warning turned off in fitcameramodel.m\n')
% warning('off','optim:lsqncommon:AlgorithmConflict')

% Variable Initialization
a_cam1=[];
a_cam2=[];
aXcam1=[];
aXcam2=[];
aYcam1=[];
aYcam2=[];

% This is a switch that will allow for compulation of the camera pinhole
% coeff. (12) at the same time as the other coeff..  Otherwise it split it
% into 3 steps.
% 1 -> for simulatanious compulation 
% 0 -> for split
allsametime=0;

% Initalization of coeff. for the different models.
if modeltype==1
    a=[100; ones(15,1)];        % initial guesss for solver
elseif modeltype==2
    a=[100; ones(18,1)];        % initial guesss for solver
else    % camera pinhole model
    
end

% Defines the percision of the outputs in the fit table in the GUI.
printformat='%5.4f';        % for converge box message
printformat1='%7.6f';

% lsqnonlin outputs
% x,resnorm,residual,exitflag,output,lambda,jacobian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If either of the cubic XY camera models are going to be used.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ((modeltype==1) || (modeltype==2))
    
    % Set model type
    alldata.orderz=modeltype;      
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get coefficients for X, cam1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    alldata.allxdata=allx1data;      
    alldata.allXdata=allX1data(:,1);
    % fit X for cam. 1 to all x data for cam1
    [aXcam1,resnorm,f_resid,exitflag,~,~,jac1]=lsqnonlin(@(a) poly_3xy_123z(a,alldata),...
        a,[],[],optionsls);   %#ok<ASGLU>
    rmsX1=sqrt(mean((f_resid).^2));%sqrt(resnorm/length(alldata.allXdata));
    corrcoef=1-sqrt(mean((f_resid./alldata.allXdata).^2));
    if ((exitflag ~= -1) && (exitflag ~= -2))
        msgX1=[' yes (' num2str(exitflag) ')     ' num2str(rmsX1,printformat) '       ' num2str(corrcoef,printformat1)];
    else
        msgX1=[' no (' num2str(exitflag) ')     ' num2str(rmsX1,printformat) '       ' num2str(corrcoef,printformat1)];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get coefficients for Y, cam1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    alldata.allxdata=allx1data;      
    alldata.allXdata=allX1data(:,2);
    % fit Y for cam. 1 to all x data for cam1
    [aYcam1,resnorm,f_resid,exitflag,~,~,jac2]=lsqnonlin(@(a) poly_3xy_123z(a,alldata),...
        a,[],[],optionsls);     %#ok<ASGLU>
    rmsY1=sqrt(mean((f_resid).^2));%sqrt(resnorm/length(alldata.allXdata));
    corrcoef=1-sqrt(mean((f_resid./alldata.allXdata).^2));
    if ((exitflag ~= -1) && (exitflag ~= -2))
        msgY1=[' yes (' num2str(exitflag) ')     ' num2str(rmsY1,printformat) '       ' num2str(corrcoef,printformat1)];
    else
        msgY1=[' no (' num2str(exitflag) ')     ' num2str(rmsY1,printformat) '       ' num2str(corrcoef,printformat1)];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get coefficients for X, cam2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    alldata.allxdata=allx2data;      
    alldata.allXdata=allX2data(:,1);
    % fit X for cam. 2 to all x data for cam1
    [aXcam2,resnorm,f_resid,exitflag,~,~,jac3]=lsqnonlin(@(a) poly_3xy_123z(a,alldata),...
        a,[],[],optionsls);   %#ok<ASGLU>  
    rmsX2=sqrt(mean((f_resid).^2));%sqrt(resnorm/length(alldata.allXdata));
    corrcoef=1-sqrt(mean((f_resid./alldata.allXdata).^2));
    if ((exitflag ~= -1) && (exitflag ~= -2))
        msgX2=[' yes (' num2str(exitflag) ')     ' num2str(rmsX2,printformat) '       ' num2str(corrcoef,printformat1)];
    else
        msgX2=[' no (' num2str(exitflag) ')     ' num2str(rmsX2,printformat) '       ' num2str(corrcoef,printformat1)];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get coefficients for Y, cam2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    alldata.allxdata=allx2data;      
    alldata.allXdata=allX2data(:,2);
    % fit Y for cam. 2 to all x data for cam1
    [aYcam2,resnorm,f_resid,exitflag,~,~,jac4]=lsqnonlin(@(a) poly_3xy_123z(a,alldata),...
        a,[],[],optionsls);   %#ok<ASGLU>  
    rmsY2=sqrt(mean((f_resid).^2));%sqrt(resnorm/length(alldata.allXdata));
    corrcoef=1-sqrt(mean((f_resid./alldata.allXdata).^2));
    if ((exitflag ~= -1) && (exitflag ~= -2))
        msgY2=[' yes (' num2str(exitflag) ')     ' num2str(rmsY2,printformat) '       ' num2str(corrcoef,printformat1)];
    else
        msgY2=[' no (' num2str(exitflag) ')     ' num2str(rmsY2,printformat) '       ' num2str(corrcoef,printformat1)];
    end
    
%     keyboard;
    
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % If modeltype =3 (camera pinhole, direct linear transformation)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if allsametime         % this way calculates a1 thru a12 all at once
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Determination of Polynomial Coef. for camera 1 using the Pinhole
        % camera model
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        a=500*ones(11,1);
        alldata.allxdata=allx1data;
        alldata.allXdata=allX1data;
        [a_cam1,resnorm,f_resid,exitflag]=lsqnonlin(@(a) campinhole_LSfun(a,alldata),...
            a,[],[],optionsls);  % fit X for cam. 1 to all x data for cam1
        %[rows cols]=size(alldata.allxdata);
        %rmsX=sqrt(resnorm/rows/2);
        %corrcoef=1-sqrt(mean((f_resid./(reshape(alldata.allXdata,rows*2,1))).^2));
        %corrcoef=1-sqrt(mean((f_resid./([reshape(alldata.allXdata,rows*2,1)])).^2));
        Xd=alldata.allXdata(:,1);
        Yd=alldata.allXdata(:,2);
        xd=alldata.allxdata(:,1);
        yd=alldata.allxdata(:,2);
        zd=alldata.allxdata(:,3);
        ad=a_cam1;
        res1=(ad(1).*xd+ad(2).*yd+ad(3).*zd+ad(4))./(ad(9).*xd+ad(10).*yd+ad(11).*zd+1)-Xd;
        res2=(ad(5).*xd+ad(6).*yd+ad(7).*zd+ad(8))./(ad(9).*xd+ad(10).*yd+ad(11).*zd+1)-Yd;
        rmsX=sqrt(mean([res1;res2].^2));
        corrcoef=1-sqrt(mean([res1./Xd;res2./Yd].^2));
        if ((exitflag ~= -1) && (exitflag ~= -2))
            msgX1=[' yes (' num2str(exitflag) ')     ' num2str(rmsX,printformat) '       ' num2str(corrcoef,printformat1)];
            msgY1=[' yes (' num2str(exitflag) ')     ' num2str(rmsX,printformat) '       ' num2str(corrcoef,printformat1)];
        else
            msgX1=[' no (' num2str(exitflag) ')     ' num2str(rmsX,printformat) '       ' num2str(corrcoef,printformat1)];
            msgY1=[' no (' num2str(exitflag) ')     ' num2str(rmsX,printformat) '       ' num2str(corrcoef,printformat1)];
        end
        a_cam1=[a_cam1;1];      % add in a34 (or a12).  this was constrained to be 1 in the solver
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Determination of Polynomial Coef. for camera 2 using the Pinhole
        % camera model
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         a=500*ones(11,1);
        alldata.allxdata=allx2data;
        alldata.allXdata=allX2data;
        [a_cam2,resnorm,f_resid,exitflag]=lsqnonlin(@(a) campinhole_LSfun(a,alldata),...
            a,[],[],optionsls);  % fit X for cam. 1 to all x data for cam1
        %[rows cols]=size(alldata.allxdata);
        %rmsX=sqrt(resnorm/rows/2);
        %corrcoef=1-sqrt(mean((f_resid./(reshape(alldata.allXdata,rows*2,1))).^2));
        %corrcoef=1-sqrt(mean((f_resid./([reshape(alldata.allXdata,rows*2,1)])).^2));
        Xd=alldata.allXdata(:,1);
        Yd=alldata.allXdata(:,2);
        xd=alldata.allxdata(:,1);
        yd=alldata.allxdata(:,2);
        zd=alldata.allxdata(:,3);
        ad=a_cam2;
        res1=(ad(1).*xd+ad(2).*yd+ad(3).*zd+ad(4))./(ad(9).*xd+ad(10).*yd+ad(11).*zd+1)-Xd;
        res2=(ad(5).*xd+ad(6).*yd+ad(7).*zd+ad(8))./(ad(9).*xd+ad(10).*yd+ad(11).*zd+1)-Yd;
        rmsX=sqrt(mean([res1;res2].^2));
        corrcoef=1-sqrt(mean([res1./Xd;res2./Yd].^2));
        if ((exitflag ~= -1) && (exitflag ~= -2))
            msgX2=[' yes (' num2str(exitflag) ')     ' num2str(rmsX,printformat) '       ' num2str(corrcoef,printformat1)];
            msgY2=[' yes (' num2str(exitflag) ')     ' num2str(rmsX,printformat) '       ' num2str(corrcoef,printformat1)];
        else
            msgX2=[' no (' num2str(exitflag) ')     ' num2str(rmsX,printformat) '       ' num2str(corrcoef,printformat1)];
            msgY2=[' no (' num2str(exitflag) ')     ' num2str(rmsX,printformat) '       ' num2str(corrcoef,printformat1)];
        end
        a_cam2=[a_cam2;1];
        
        
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This way calculates a1 thru a12 on three separate times, they 
        % turn out to be the same
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        a=ones(1,4);
        alldata.allxdata=allx1data;
        alldata.allXdata=allX1data(:,1);
        % fit X for cam. 1 to all x data for cam1
        [a1_cam1,resnorm1,f_resid1,~]=lsqnonlin(@(a) campinholemod_LSfun(a,alldata),...
            a,[],[],optionsls);  
        
        a=ones(1,4);
        alldata.allxdata=allx1data;
        alldata.allXdata=allX1data(:,2);
        % fit X for cam. 1 to all y data for cam1
        [a2_cam1,resnorm2,f_resid2,~]=lsqnonlin(@(a) campinholemod_LSfun(a,alldata),...
            a,[],[],optionsls);  
        
        a=ones(1,4);
        alldata.allxdata=allx1data;
        [rows,~]=size(alldata.allxdata);
        alldata.allXdata=ones(rows,1);
        % fit X for cam. 1 to all z data for cam1
        [a3_cam1,resnorm3,f_resid3,exitflag]=lsqnonlin(@(a) campinholemod_LSfun(a,alldata),...
            a,[],[],optionsls);  
        
        rmsX=sqrt((resnorm1+resnorm2+resnorm3)/3/rows);
        corrcoef=1-sqrt(mean([(f_resid1./(allX1data(:,1))).^2;(f_resid2./(allX1data(:,2))).^2;(f_resid3./(ones(rows,1))).^2]));
        
        msgX1=[' yes (' num2str(exitflag) ')     ' num2str(rmsX,printformat) '       ' num2str(corrcoef,printformat1)];
        msgY1=[' yes (' num2str(exitflag) ')     ' num2str(rmsX,printformat) '       ' num2str(corrcoef,printformat1)];
        
        %[f_resid1 f_resid2 f_resid3]
        
        a=ones(1,4);
        alldata.allxdata=allx2data;
        alldata.allXdata=allX2data(:,1);
        % fit X for cam. 2 to all x data for cam1
        [a1_cam2,resnorm1,f_resid1,~]=lsqnonlin(@(a) campinholemod_LSfun(a,alldata),...
            a,[],[],optionsls); 
        
        a=ones(1,4);
        alldata.allxdata=allx2data;
        alldata.allXdata=allX2data(:,2);
        % fit X for cam. 2 to all y data for cam1
        [a2_cam2,resnorm2,f_resid2,~]=lsqnonlin(@(a) campinholemod_LSfun(a,alldata),...
            a,[],[],optionsls);
        
        a=ones(1,4);
        alldata.allxdata=allx2data;
        [rows,~]=size(alldata.allxdata);
        alldata.allXdata=ones(rows,1);
        % fit X for cam. 2 to all z data for cam1
        [a3_cam2,resnorm3,f_resid3,exitflag]=lsqnonlin(@(a) campinholemod_LSfun(a,alldata),...
            a,[],[],optionsls);  
        
        rmsX=sqrt((resnorm1+resnorm2+resnorm3)/3/rows);
        corrcoef=1-sqrt(mean([(f_resid1./(allX2data(:,1))).^2;(f_resid2./(allX2data(:,2))).^2;(f_resid3./(ones(rows,1))).^2]));
        
        msgX2=[' yes (' num2str(exitflag) ')     ' num2str(rmsX,printformat) '       ' num2str(corrcoef,printformat1)];
        msgY2=[' yes (' num2str(exitflag) ')     ' num2str(rmsX,printformat) '       ' num2str(corrcoef,printformat1)];
        
        % [f_resid1 f_resid2 f_resid3]
        
        % Combine the camera coef for each camera
        a_cam1=[a1_cam1; a2_cam1; a3_cam1];
        a_cam2=[a1_cam2; a2_cam2; a3_cam2];
    end
end
convergemessage={msgX1;msgY1;msgX2;msgY2};
% Approximate camera angles calculated for each calibration Mapping
% function using the ratio of the coefficients.
% if modeltype==1 || modeltype==2
%     
%     [aXcam1 aYcam1 aXcam2 aYcam2]
%     fprintf('\n Approximate Camera Angles: \n Alpha1  Alpha2  Beta1 Beta2\n');
%     alpha1=atand((aYcam1(4)*aXcam1(3) - aYcam1(3)*aXcam1(4))/(aYcam1(3)*aXcam1(2) - aYcam1(2)*aXcam1(3)));
%     alpha2=atand((aYcam2(4)*aXcam2(3) - aYcam2(3)*aXcam2(4))/(aYcam2(3)*aXcam2(2) - aYcam2(2)*aXcam2(3)));
%     beta1=atand((aYcam1(4)*aXcam1(2) - aYcam1(2)*aXcam1(4))/(aYcam1(2)*aXcam1(3) - aYcam1(3)*aXcam1(2)));
%     beta2=atand((aYcam2(4)*aXcam2(2) - aYcam2(2)*aXcam2(4))/(aYcam2(2)*aXcam2(3) - aYcam2(3)*aXcam2(2)));
%     [alpha1 alpha2 beta1 beta2]
%     
% elseif modeltype==3
%     [a_cam1;a_cam2]
%     aXcam1=[a_cam1(1,4) a_cam1(1,1) a_cam1(1,2) a_cam1(1,3)]';
%     aYcam1=[a_cam1(2,4) a_cam1(2,1) a_cam1(2,2) a_cam1(2,3)]';
%     aXcam2=[a_cam2(1,4) a_cam2(1,1) a_cam2(1,2) a_cam2(1,3)]';
%     aYcam2=[a_cam2(2,4) a_cam2(2,1) a_cam2(2,2) a_cam2(2,3)]';
%     [aXcam1 aYcam1 aXcam2 aYcam2]
%     
%     fprintf('\n Approximate Camera Angles: \n Alpha1  Alpha2  Beta1 Beta2\n');
%     alpha1=atand((aYcam1(4)*aXcam1(3) - aYcam1(3)*aXcam1(4))/(aYcam1(3)*aXcam1(2) - aYcam1(2)*aXcam1(3)));
%     alpha2=atand((aYcam2(4)*aXcam2(3) - aYcam2(3)*aXcam2(4))/(aYcam2(3)*aXcam2(2) - aYcam2(2)*aXcam2(3)));
%     beta1=atand((aYcam1(4)*aXcam1(2) - aYcam1(2)*aXcam1(4))/(aYcam1(2)*aXcam1(3) - aYcam1(3)*aXcam1(2)));
%     beta2=atand((aYcam2(4)*aXcam2(2) - aYcam2(2)*aXcam2(4))/(aYcam2(2)*aXcam2(3) - aYcam2(3)*aXcam2(2)));
%     [alpha1 alpha2 beta1 beta2]
%     
%     aXcam1=[];aYcam1=[];aXcam2=[];aYcam2=[];
%     
% end
% save('aXcam1.mat','aXcam1');
% save('aXcam2.mat','aXcam2');
% save('aYcam1.mat','aYcam1');
% save('aYcam2.mat','aYcam2');
% Turning the warning back on
%warning('on','optim:lsqncommon:AlgorithmConflict')

% Calculating fit coefficient Uncertainty
Np=size(allx1data,1);
Na=length(aXcam1);
Ndof=Np-Na;

SigresX1=sqrt(Np/Ndof)*rmsX1*eye(Na);
SigresY1=sqrt(Np/Ndof)*rmsY1*eye(Na);
SigresX2=sqrt(Np/Ndof)*rmsX2*eye(Na);
SigresY2=sqrt(Np/Ndof)*rmsY2*eye(Na);

% Cov1=inv(jac1'*jac1);
% Cov2=inv(jac2'*jac2);
% Cov3=inv(jac3'*jac3);
% Cov4=inv(jac4'*jac4);

Cov1=(jac1'*jac1);
Cov2=(jac2'*jac2);
Cov3=(jac3'*jac3);
Cov4=(jac4'*jac4);

Uxcam1=sqrt(diag(Cov1\SigresX1.^2));
Uycam1=sqrt(diag(Cov2\SigresY1.^2));
Uxcam2=sqrt(diag(Cov3\SigresX2.^2));
Uycam2=sqrt(diag(Cov4\SigresY2.^2));

Uncal=[Uxcam1 Uycam1 Uxcam2 Uycam2];

end
