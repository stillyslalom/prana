function [Un_u,Un_v,Un_w,JU,JV,JW]=stereo_uncertainty_propagation(Unu1,Unv1,Unu2,Unv2,U1,V1,U2,V2,W,Un_alpha1,Un_alpha2,Un_beta1,Un_beta2,tanalpha1,tanalpha2,tanbeta1,tanbeta2,mx,my)
% This function calculates the uncertainty in the stereo velocity
% components
%
% Input Variables
% Nx X Ny grid point array on which individual camera velocity fields are
% evaluated
% Unu1= Nx X Ny array of uncertainty in camera 1 U component of velocity using any planar uncertainty method like IM or CS
% Unv1= Nx X Ny array of uncertainty in camera 1 V component of velocity using any planar uncertainty method like IM or CS
% Unu2= Nx X Ny array of uncertainty in camera 2 U component of velocity using any planar uncertainty method like IM or CS
% Unv2= Nx X Ny array of uncertainty in camera 2 V component of velocity using any planar uncertainty method like IM or CS
% W= Out of plane velocity velocity component
%
% tanalpha1,tanalpha2,tanbeta1,tanbeta2= tangent of the stereo angles alpha
% (angle in x-z plane) and beta (angle in y-z plane) for camera 1 and 2
%
% Un_alpha1,Un_alpha2,Un_beta1,Un_beta2= the uncertainty in stereo
% angles alpha1, beta1, alpha2, beta2
%
% The angles and its uncertainties are calculated previously in
% stereo_angle_uncertainty.m
%
% mx,my = magnifications in mm(or physical unit)/pixel for the dewarped
% common grid in x and y direction respectively.
%
% Output Variables
% Un_u,Un_v and Un_w are the uncertainties in the stereo velocity
% components u, v and w.
%
%JU, JV and JW = Sensitivity coefficients of U, V and W uncertainty propagation equation as given in Table 1 of the stereo uncertainty paper   

%written by Sayantan Bhattacharya on 01/12/2016


%Get the stereo angles from the tangent of the angles
A1=atan(tanalpha1);
A2=atan(tanalpha2);
B1=atan(tanbeta1);
B2=atan(tanbeta2);


if (max(abs(tanalpha1(:)))+ max(abs(tanalpha2(:))))> (max(abs(tanbeta1(:)))+ max(abs(tanbeta2(:)))) % For horizontal configuration (X-Z plane)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%Uncertainty horizontal camera%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% calculating coefficients for U component
    %Sensitivity coeeficients for U-component uncertainty propagation
    Udu1= tanalpha2./(tanalpha2-tanalpha1);
    Udu2= -tanalpha1./(tanalpha2-tanalpha1);
    Uda1= ((U1-U2).*sin(A2).*cos(A2))./(sin(A2-A1)).^2;
    Uda2= ((U2-U1).*sin(A1).*cos(A1))./(sin(A2-A1)).^2;
    
    lm=zeros(size(Unu1(:)));
    SigU=zeros(4,4,size(lm,1));
    JU=zeros(1,4,size(lm,1));
    
    Un_u=zeros(size(lm));Un_v=zeros(size(lm));Un_w=zeros(size(lm));
    
    %Calculating U jacobian
    JU(1,1,:)=Udu1(:); JU(1,2,:)=Udu2(:); JU(1,3,:)=Uda1(:); JU(1,4,:)=Uda2(:);
    
    %variances or diagonal terms of covariance matrix
    su1=Unu1(:).^2;
    su2=Unu2(:).^2;
    su3=Un_alpha1(:).^2;
    su4=Un_alpha2(:).^2;
   
    %correlation coefficients for elemental uncertainties
    %This is set to zero and thus the covariance terms are neglected for
    %now
    ru12=0;ru13=0;ru14=0;ru23=0;ru24=0;ru34=0;

    %Covariance matrix for U
    SigU(1,1,:)=su1;                    SigU(1,2,:)=ru12.*sqrt(su1.*su2);   SigU(1,3,:)=ru13.*sqrt(su1.*su3);   SigU(1,4,:)=ru14.*sqrt(su1.*su4);
    SigU(2,1,:)=ru12.*sqrt(su1.*su2);   SigU(2,2,:)=su2;                    SigU(2,3,:)=ru23.*sqrt(su2.*su3);   SigU(2,4,:)=ru24.*sqrt(su2.*su4); 
    SigU(3,1,:)=ru13.*sqrt(su1.*su3);   SigU(3,2,:)=ru23.*sqrt(su2.*su3);   SigU(3,3,:)=su3;                    SigU(3,4,:)=ru34.*sqrt(su3.*su4); 
    SigU(4,1,:)=ru14.*sqrt(su1.*su4);   SigU(4,2,:)=ru24.*sqrt(su2.*su4);   SigU(4,3,:)=ru34.*sqrt(su3.*su4);   SigU(4,4,:)=su4; 

   
    %% calculating coefficients for W component
    %Sensitivity coeeficients for W-component uncertainty propagation
    Wdu1= 1./(tanalpha2-tanalpha1);
    Wdu2= -1./(tanalpha2-tanalpha1);
    Wda1= ((U1-U2).*cos(A2).^2)./(sin(A2-A1)).^2;
    Wda2= -((U1-U2).*cos(A1).^2)./(sin(A2-A1)).^2;
    
    SigW=zeros(4,4,size(lm,1));
    JW=zeros(1,4,size(lm,1));
    
    %Calculating W jacobian
    JW(1,1,:)=Wdu1(:); JW(1,2,:)=Wdu2(:); JW(1,3,:)=Wda1(:); JW(1,4,:)=Wda2(:);
    
    %%variances or diagonal terms of covariance matrix
    sw1=Unu1(:).^2;
    sw2=Unu2(:).^2;
    sw3=Un_alpha1(:).^2;
    sw4=Un_alpha2(:).^2;
    
    %correlation coefficients for elemental uncertainties
    %This is set to zero and thus the covariance terms are neglected for
    %now
    rw12=0;rw13=0;rw14=0;rw23=0;rw24=0;rw34=0;
    
    %Covariance matrix for W
    SigW(1,1,:)=sw1;                    SigW(1,2,:)=rw12.*sqrt(sw1.*sw2);   SigW(1,3,:)=rw13.*sqrt(sw1.*sw3);   SigW(1,4,:)=rw14.*sqrt(sw1.*sw4);
    SigW(2,1,:)=rw12.*sqrt(sw1.*sw2);   SigW(2,2,:)=sw2;                    SigW(2,3,:)=rw23.*sqrt(sw2.*sw3);   SigW(2,4,:)=rw24.*sqrt(sw2.*sw4); 
    SigW(3,1,:)=rw13.*sqrt(sw1.*sw3);   SigW(3,2,:)=rw23.*sqrt(sw2.*sw3);   SigW(3,3,:)=sw3;                    SigW(3,4,:)=rw34.*sqrt(sw3.*sw4); 
    SigW(4,1,:)=rw14.*sqrt(sw1.*sw4);   SigW(4,2,:)=rw24.*sqrt(sw2.*sw4);   SigW(4,3,:)=rw34.*sqrt(sw3.*sw4);   SigW(4,4,:)=sw4; 
    
    for i=1:length(lm)
        Un_u(i)=sqrt(JU(:,:,i)*SigU(:,:,i)*JU(:,:,i)');
        Un_w(i)=sqrt(JW(:,:,i)*SigW(:,:,i)*JW(:,:,i)');
    end
%     keyboard;
    %% calculating coefficients for V component
    %Sensitivity coeeficients for V-component uncertainty propagation
    Vdv1= 0.5;
    Vdv2= 0.5;
    Vdb1= (mx/my)*(W./2).*(1./(cos(B1).^2));
    Vdb2= (mx/my)*(W./2).*(1./(cos(B2).^2));
    VdW = (mx/my)*0.5.*(tanbeta1+tanbeta2);
    
    
    SigV=zeros(5,5,size(lm,1));
    JV=zeros(1,5,size(lm,1));
    
    %Calculating W jacobian
    JV(1,1,:)=Vdv1(:); JV(1,2,:)=Vdv2(:); JV(1,3,:)=Vdb1(:); JV(1,4,:)=Vdb2(:);  JV(1,5,:)=VdW(:);
    
    %%variances or diagonal terms of covariance matrix
    sv1=Unv1(:).^2;
    sv2=Unv2(:).^2;
    sv3=Un_beta1(:).^2;
    sv4=Un_beta2(:).^2;
    sv5=Un_w(:).^2;
    
    %correlation coefficients for elemental uncertainties
    %This is set to zero and thus the covariance terms are neglected for
    %now
    rv12=0;rv13=0;rv14=0;rv15=0;rv23=0;rv24=0;rv25=0;rv34=0;rv35=0;rv45=0;

    %Covariance matrix for V
    SigV(1,1,:)=sv1;                    SigV(1,2,:)=rv12.*sqrt(sv1.*sv2);   SigV(1,3,:)=rv13.*sqrt(sv1.*sv3);   SigV(1,4,:)=rv14.*sqrt(sv1.*sv4);   SigV(1,5,:)=rv15.*sqrt(sv1.*sv5);
    SigV(2,1,:)=rv12.*sqrt(sv1.*sv2);   SigV(2,2,:)=sv2;                    SigV(2,3,:)=rv23.*sqrt(sv2.*sv3);   SigV(2,4,:)=rv24.*sqrt(sv2.*sv4);   SigV(2,5,:)=rv25.*sqrt(sv2.*sv5); 
    SigV(3,1,:)=rv13.*sqrt(sv1.*sv3);   SigV(3,2,:)=rv23.*sqrt(sv2.*sv3);   SigV(3,3,:)=sv3;                    SigV(3,4,:)=rv34.*sqrt(sv3.*sv4);   SigV(3,5,:)=rv35.*sqrt(sv3.*sv5); 
    SigV(4,1,:)=rv14.*sqrt(sv1.*sv4);   SigV(4,2,:)=rv24.*sqrt(sv2.*sv4);   SigV(4,3,:)=rv34.*sqrt(sv3.*sv4);   SigV(4,4,:)=sv4;                    SigV(4,5,:)=rv45.*sqrt(sv4.*sv5);
    SigV(5,1,:)=rv15.*sqrt(sv1.*sv5);   SigV(5,2,:)=rv25.*sqrt(sv2.*sv5);   SigV(5,3,:)=rv35.*sqrt(sv3.*sv5);   SigV(5,4,:)=rv45.*sqrt(sv4.*sv5);   SigV(5,5,:)=sv5;
    
    % The uncertainty propagation equation
    for i=1:length(lm)
        Un_v(i)=sqrt(JV(:,:,i)*SigV(:,:,i)*JV(:,:,i)');
    end
    
    Un_u=reshape(Un_u,size(Unu1,1),size(Unu1,2));
    Un_v=reshape(Un_v,size(Unu1,1),size(Unu1,2));
    Un_w=reshape(Un_w,size(Unu1,1),size(Unu1,2));

    
    
elseif (max(abs(tanalpha1(:)))+ max(abs(tanalpha2(:))))< (max(abs(tanbeta1(:)))+ max(abs(tanbeta2(:)))) % For Vertical Configuration (Y-Z plane)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%Uncertainty vertical camera%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% calculating coefficients for V component
    %Sensitivity coeeficients for V-component uncertainty propagation
    Vdv1= tanbeta2./(tanbeta2-tanbeta1);
    Vdv2= -tanbeta1./(tanbeta2-tanbeta1);
    Vdb1= ((V1-V2).*sin(B2).*cos(B2))./(sin(B2-B1)).^2;
    Vdb2= ((V2-V1).*sin(B1).*cos(B1))./(sin(B2-B1)).^2;
    
    lm=zeros(size(Unv1(:)));
    SigV=zeros(4,4,size(lm,1));
    JV=zeros(1,4,size(lm,1));
    
    Un_u=zeros(size(lm));Un_v=zeros(size(lm));Un_w=zeros(size(lm));
    
    %Calculating V jacobian
    JV(1,1,:)=Vdv1(:); JV(1,2,:)=Vdv2(:); JV(1,3,:)=Vdb1(:); JV(1,4,:)=Vdb2(:);
    %keyboard;
    
    %variances or diagonal terms of covariance matrix
    sv1=Unv1(:).^2;
    sv2=Unv2(:).^2;
    sv3=Un_beta1(:).^2;
    sv4=Un_beta2(:).^2;
    
    %correlation coefficients for elemental uncertainties
    %This is set to zero and thus the covariance terms are neglected for
    %now
    rv12=0;rv13=0;rv14=0;rv23=0;rv24=0;rv34=0;
    
    %Covariance matrix for V
    SigV(1,1,:)=sv1;                    SigV(1,2,:)=rv12.*sqrt(sv1.*sv2);   SigV(1,3,:)=rv13.*sqrt(sv1.*sv3);   SigV(1,4,:)=rv14.*sqrt(sv1.*sv4);
    SigV(2,1,:)=rv12.*sqrt(sv1.*sv2);   SigV(2,2,:)=sv2;                    SigV(2,3,:)=rv23.*sqrt(sv2.*sv3);   SigV(2,4,:)=rv24.*sqrt(sv2.*sv4); 
    SigV(3,1,:)=rv13.*sqrt(sv1.*sv3);   SigV(3,2,:)=rv23.*sqrt(sv2.*sv3);   SigV(3,3,:)=sv3;                    SigV(3,4,:)=rv34.*sqrt(sv3.*sv4); 
    SigV(4,1,:)=rv14.*sqrt(sv1.*sv4);   SigV(4,2,:)=rv24.*sqrt(sv2.*sv4);   SigV(4,3,:)=rv34.*sqrt(sv3.*sv4);   SigV(4,4,:)=sv4; 
    

    
    %% calculating coefficients for W component
    %Sensitivity coeeficients for W-component uncertainty propagation
    Wdv1= 1./(tanbeta2-tanbeta1);
    Wdv2= -1./(tanbeta2-tanbeta1);
    Wdb1= ((V1-V2).*cos(B2).^2)./(sin(B2-B1)).^2;
    Wdb2= -((V1-V2).*cos(B1).^2)./(sin(B2-B1)).^2;
    
    SigW=zeros(4,4,size(lm,1));
    JW=zeros(1,4,size(lm,1));
    
    %Calculating W jacobian
    JW(1,1,:)=Wdv1(:); JW(1,2,:)=Wdv2(:); JW(1,3,:)=Wdb1(:); JW(1,4,:)=Wdb2(:);
    
    %%variances or diagonal terms of covariance matrix
    sw1=Unv1(:).^2;
    sw2=Unv2(:).^2;
    sw3=Un_beta1(:).^2;
    sw4=Un_beta2(:).^2;
    
    %correlation coefficients for elemental uncertainties
    %This is set to zero and thus the covariance terms are neglected for
    %now
    rw12=0;rw13=0;rw14=0;rw23=0;rw24=0;rw34=0;

    
    %Covariance matrix for W
    SigW(1,1,:)=sw1;                    SigW(1,2,:)=rw12.*sqrt(sw1.*sw2);   SigW(1,3,:)=rw13.*sqrt(sw1.*sw3);   SigW(1,4,:)=rw14.*sqrt(sw1.*sw4);
    SigW(2,1,:)=rw12.*sqrt(sw1.*sw2);   SigW(2,2,:)=sw2;                    SigW(2,3,:)=rw23.*sqrt(sw2.*sw3);   SigW(2,4,:)=rw24.*sqrt(sw2.*sw4); 
    SigW(3,1,:)=rw13.*sqrt(sw1.*sw3);   SigW(3,2,:)=rw23.*sqrt(sw2.*sw3);   SigW(3,3,:)=sw3;                    SigW(3,4,:)=rw34.*sqrt(sw3.*sw4); 
    SigW(4,1,:)=rw14.*sqrt(sw1.*sw4);   SigW(4,2,:)=rw24.*sqrt(sw2.*sw4);   SigW(4,3,:)=rw34.*sqrt(sw3.*sw4);   SigW(4,4,:)=sw4; 
    
    for i=1:length(lm)
        Un_v(i)=sqrt(JV(:,:,i)*SigV(:,:,i)*JV(:,:,i)');
        Un_w(i)=sqrt(JW(:,:,i)*SigW(:,:,i)*JW(:,:,i)');
    end
    %keyboard;
    %% calculating coefficients for U component
    %Sensitivity coeeficients for U-component uncertainty propagation
    Udu1= 0.5;
    Udu2= 0.5;
    Uda1= (my/mx)*(W./2).*(1./(cos(A1).^2));
    Uda2= (my/mx)*(W./2).*(1./(cos(A2).^2));
    UdW = (my/mx)*0.5.*(tanalpha1+tanalpha2);
    
    
    SigU=zeros(5,5,size(lm,1));
    JU=zeros(1,5,size(lm,1));
    
    %Calculating U jacobian
    JU(1,1,:)=Udu1(:); JU(1,2,:)=Udu2(:); JU(1,3,:)=Uda1(:); JU(1,4,:)=Uda2(:);  JU(1,5,:)=UdW(:);
    
    %%variances or diagonal terms of covariance matrix
    su1=Unu1(:).^2;
    su2=Unu2(:).^2;
    su3=Un_alpha1(:).^2;
    su4=Un_alpha2(:).^2;
    su5=Un_w(:).^2;
    
    %correlation coefficients for elemental uncertainties
    %This is set to zero and thus the covariance terms are neglected for
    %now
    ru12=0;ru13=0;ru14=0;ru15=0;ru23=0;ru24=0;ru25=0;ru34=0;ru35=0;ru45=0;
    
    %Couariance matrix for U
    SigU(1,1,:)=su1;                    SigU(1,2,:)=ru12.*sqrt(su1.*su2);   SigU(1,3,:)=ru13.*sqrt(su1.*su3);   SigU(1,4,:)=ru14.*sqrt(su1.*su4);   SigU(1,5,:)=ru15.*sqrt(su1.*su5);
    SigU(2,1,:)=ru12.*sqrt(su1.*su2);   SigU(2,2,:)=su2;                    SigU(2,3,:)=ru23.*sqrt(su2.*su3);   SigU(2,4,:)=ru24.*sqrt(su2.*su4);   SigU(2,5,:)=ru25.*sqrt(su2.*su5); 
    SigU(3,1,:)=ru13.*sqrt(su1.*su3);   SigU(3,2,:)=ru23.*sqrt(su2.*su3);   SigU(3,3,:)=su3;                    SigU(3,4,:)=ru34.*sqrt(su3.*su4);   SigU(3,5,:)=ru35.*sqrt(su3.*su5); 
    SigU(4,1,:)=ru14.*sqrt(su1.*su4);   SigU(4,2,:)=ru24.*sqrt(su2.*su4);   SigU(4,3,:)=ru34.*sqrt(su3.*su4);   SigU(4,4,:)=su4;                    SigU(4,5,:)=ru45.*sqrt(su4.*su5);
    SigU(5,1,:)=ru15.*sqrt(su1.*su5);   SigU(5,2,:)=ru25.*sqrt(su2.*su5);   SigU(5,3,:)=ru35.*sqrt(su3.*su5);   SigU(5,4,:)=ru45.*sqrt(su4.*su5);   SigU(5,5,:)=su5;
    
    % The uncertainty propagation equation
    for i=1:length(lm)
        Un_u(i)=sqrt(JU(:,:,i)*SigU(:,:,i)*JU(:,:,i)');
    end
    
    Un_u=reshape(Un_u,size(Unu1,1),size(Unu1,2));
    Un_v=reshape(Un_v,size(Unu1,1),size(Unu1,2));
    Un_w=reshape(Un_w,size(Unu1,1),size(Unu1,2));
    

end


end