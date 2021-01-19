function [Un_alpha1,Un_alpha2,Un_beta1,Un_beta2]=stereo_uncertainty_in_map_func_grad(calmat,xg,yg,zg,Uncalcoeff,unwx,unwy,unwz,modeltype)
% This function calculates the uncertainty in the mapping function gradients and subsequesntly the angle uncertainty
%written by Sayantan Bhattacharya January 2016    
[r,c]=size(xg);
dFdx1=zeros(r,c,4);       % the 3rd dimention corresponds to dFdx1 for (X1,Y1,X2,Y2)
dFdx2=zeros(r,c,4);
dFdx3=zeros(r,c,4);
Un_dFdx1=zeros(r,c,4);       % the 3rd dimention corresponds to dFdx1 for (X1,Y1,X2,Y2)
Un_dFdx2=zeros(r,c,4);
Un_dFdx3=zeros(r,c,4);

% Mapping function gradient
[~,~,dFdx1(:,:,1:2),dFdx2(:,:,1:2),dFdx3(:,:,1:2)]=calculate_stereo_angle(calmat(:,1:2),xg,yg,zg,modeltype);
[~,~,dFdx1(:,:,3:4),dFdx2(:,:,3:4),dFdx3(:,:,3:4)]=calculate_stereo_angle(calmat(:,3:4),xg,yg,zg,modeltype);


for gg=1:4
    if modeltype==2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Mapping the camera coord. to the World Coord. using 1st order z
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        a=calmat(:,gg); %original mapping function coefficients
        %{    
        % Mapping function gradient
        dFdx1(:,:,gg) = a(2) + 2*a(5)*xg + a(6)*yg + a(8)*zg + 3*a(10)*xg.^2 + ...
            2*a(11)*xg.*yg + a(12)*yg.^2 + 2*a(14)*xg.*zg + a(15)*yg.*zg;
        
        dFdx2(:,:,gg) = a(3) + a(6)*xg + 2*a(7)*yg + a(9)*zg + a(11)*xg.^2 + ...
            2*a(12)*xg.*yg + 3*a(13)*yg.^2 + a(15)*xg.*zg + 2*a(16)*yg.*zg;
        
        dFdx3(:,:,gg) = a(4) + a(8)*xg + a(9)*yg + a(14)*xg.^2 + a(15)*xg.*yg + a(16)*yg.^2;
        %}

        % Second order derivative
        cfxx=2*a(5) +6*a(10)*xg + 2*a(11)*yg + 2*a(14)*zg;
        cfxy=a(6) + 2*a(11)*xg +  2*a(12)*yg + a(15)*zg;
        cfxz=a(8)+ 2*a(14)*xg + a(15)*yg;
    
        cfyx=a(6) + 2*a(11)*xg +  2*a(12)*yg + a(15)*zg;
        cfyy=2*a(7) +2*a(12)*xg + 6*a(13)*yg +  2*a(16)*zg;
        cfyz= a(9)+ a(15)*xg + 2*a(16)*yg;
    
        cfzx=a(8)+ 2*a(14)*xg + a(15)*yg;
        cfzy=a(9)+ a(15)*xg + 2*a(16)*yg;
        cfzz=0*zg;
        
        
        % Uncertainty due to cal coefficients
        a=Uncalcoeff(:,gg); % Uncertainty in mapping function coefficients
        
        T2dFdx1=sqrt((a(2)).^2 + (2*a(5)*xg).^2 + (a(6)*yg).^2 + (a(8)*zg).^2 + (3*a(10)*xg.^2).^2 + ...
            (2*a(11)*xg.*yg).^2 + (a(12)*yg.^2).^2 + (2*a(14)*xg.*zg).^2 + (a(15)*yg.*zg).^2);
    
        T2dFdx2=sqrt((a(3)).^2 + (a(6)*xg).^2 + (2*a(7)*yg).^2 + (a(9)*zg).^2 + (a(11)*xg.^2).^2 + ...
            (2*a(12)*xg.*yg).^2 + (3*a(13)*yg.^2).^2 + (a(15)*xg.*zg).^2 + (2*a(16)*yg.*zg).^2);
    
        T2dFdx3=sqrt((a(4)).^2 + (a(8)*xg).^2 + (a(9)*yg).^2 + (a(14)*xg.^2).^2 + ...
            (a(15)*xg.*yg).^2 + (a(16)*yg.^2).^2);
        
    elseif modeltype==2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Mapping the camera coord. to the World Coord. using 2nd order z
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        a=calmat(:,gg);
        %        Uc=uncoeff(:,gg);
    
        %{    
        % Mapping function gradient
        dFdx1(:,:,gg) = a(2) + 2*a(5)*xg + a(6)*yg + a(8)*zg + 3*a(11)*xg.^2 + 2*a(12)*xg.*yg + ...
            a(13)*yg.^2 + 2*a(15)*xg.*zg + a(16)*yg.*zg + a(18)*zg.^2;
    
        dFdx2(:,:,gg) = a(3) + a(6)*xg + 2*a(7)*yg + a(9)*zg + a(12)*xg.^2 + 2*a(13)*xg.*yg + ...
            3*a(14)*yg.^2 + a(16)*xg.*zg + 2*a(17)*yg.*zg + a(19)*zg.^2;
    
        dFdx3(:,:,gg) = a(4) + a(8)*xg + a(9)*yg + 2*a(10)*zg + a(15)*xg.^2 + a(16)*xg.*yg + ...
            a(17)*yg.^2 + 2*a(18)*xg.*zg + 2*a(19)*yg.*zg;
        %}
    
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
    elseif modeltype==4
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Mapping the camera coord. to the World Coord. using linear interp between cubic xy planes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        a=calmat(:,gg); %original mapping function coefficients
        %{
        dFdx1(:,:,gg) = a(2) + 2*a(5)*xg + a(6)*yg + a(8)*zg + 3*a(10)*xg.^2 + ...
            2*a(11)*xg.*yg + a(12)*yg.^2 + 2*a(14)*xg.*zg + a(15)*yg.*zg + ...
            2*a(17)*xg.*yg.*zg + a(18)*yg.^2.*zg + 3*a(19)*xg.^2.*zg;

        dFdx2(:,:,gg) = a(3) + a(6)*xg + 2*a(7)*yg + a(9)*zg + a(11)*xg.^2 + ...
            2*a(12)*xg.*yg + 3*a(13)*yg.^2 + a(15)*xg.*zg + 2*a(16)*yg.*zg + ...
            a(17)*xg.^2.*zg + 2*a(18)*xg.*yg.*zg + 3*a(20)*yg.^2.*zg;

        dFdx3(:,:,gg) = a(4) + a(8)*xg + a(9)*yg + a(14)*xg.^2 + a(15)*xg.*yg + a(16)*yg.^2 + ...
            a(17)*xg.^2.*yg + a(18)*xg.*yg.^2 + a(19)*xg.^3 + a(20)*yg.^3;
        %}
        
        % Second order derivative 
        cfxx=2*a(5) + 6*a(10)*xg + 2*a(11)*yg + 2*a(14)*zg + 2*a(17)*yg.*zg + 6*a(19)*xg.*zg;
        cfxy=a(6) + 2*a(11)*xg + 2*a(12)*yg + a(15)*zg + 2*a(17)*xg.*zg + 2*a(18)*yg.*zg;
        cfxz=a(8) + 2*a(14)*xg + a(15)*yg + 2*a(17)*xg.*yg + a(18)*yg.^2 + 3*a(19)*xg.^2;

        cfyx=a(6) + 2*a(11)*xg + 2*a(12)*yg + a(15)*zg + 2*a(17)*xg.*zg + 2*a(18)*yg.*zg;
        cfyy=2*a(7) + 2*a(12)*xg + 6*a(13)*yg + 2*a(16)*zg + 2*a(18)*xg.*zg + 6*a(20)*yg.*zg;
        cfyz=a(9) + a(15)*xg + 2*a(16)*yg + a(17)*xg.^2 + 2*a(18)*xg.*yg + 3*a(20)*yg.^2;

        cfzx=a(8) + 2*a(14)*xg + a(15)*yg + 2*a(17)*xg.*yg + a(18)*yg.^2 + 3*a(19)*xg.^2;
        cfzy=a(9) + a(15)*xg + 2*a(16)*yg + a(17)*xg.^2 + 2*a(18)*xg.*yg + 3*a(20)*yg.^2;
        cfzz=0*zg;
        
        % Uncertainty due to cal coefficients
        a=Uncalcoeff(:,gg); % Uncertainty in mapping function coefficients
        
        T2dFdx1=sqrt((a(2)).^2 + (2*a(5)*xg).^2 + (a(6)*yg).^2 + (a(8)*zg).^2 + (3*a(10)*xg.^2).^2 + ...
            (2*a(11)*xg.*yg).^2 + (a(12)*yg.^2).^2 + (2*a(14)*xg.*zg).^2 + (a(15)*yg.*zg).^2 + ...
            (2*a(17)*xg.*yg.*zg).^2 + (a(18)*yg.^2.*zg).^2 + (3*a(19)*xg.^2.*zg).^2);
    
        T2dFdx2=sqrt((a(3)).^2 + (a(6)*xg).^2 + (2*a(7)*yg).^2 + (a(9)*zg).^2 + (a(11)*xg.^2).^2 + ...
            (2*a(12)*xg.*yg).^2 + (3*a(13)*yg.^2).^2 + (a(15)*xg.*zg).^2 + (2*a(16)*yg.*zg).^2 + ...
            (a(17)*xg.^2.*zg).^2 + (2*a(18)*xg.*yg.*zg).^2 + (3*a(20)*yg.^2.*zg).^2);
    
        T2dFdx3=sqrt((a(4)).^2 + (a(8)*xg).^2 + (a(9)*yg).^2 + (a(14)*xg.^2).^2 + (a(15)*xg.*yg).^2 + (a(16)*yg.^2).^2 + ...
            (a(17)*xg.^2.*yg).^2 + (a(18)*xg.*yg.^2).^2 + (a(19)*xg.^3).^2 + (a(20)*yg.^3).^2);
    end
    
    %Uncertainty due to Ux, Uy and Uz
    T1dFdx1=sqrt((cfxx.^2).*(unwx.^2)+(cfxy.^2).*(unwy.^2)+(cfxz.^2).*(unwz.^2));
    T1dFdx2=sqrt((cfyx.^2).*(unwx.^2)+(cfyy.^2).*(unwy.^2)+(cfyz.^2).*(unwz.^2));
    T1dFdx3=sqrt((cfzx.^2).*(unwx.^2)+(cfzy.^2).*(unwy.^2)+(cfzz.^2).*(unwz.^2));
    
    %Total Uncertainty in mapping function gradient
    Un_dFdx1(:,:,gg)=sqrt(T1dFdx1.^2+T2dFdx1.^2);
    Un_dFdx2(:,:,gg)=sqrt(T1dFdx2.^2+T2dFdx2.^2);
    Un_dFdx3(:,:,gg)=sqrt(T1dFdx3.^2+T2dFdx3.^2);

    % % Uncertainty in mapping function gradient due to Ux,Uy,Uz
    % Un_dFdx1(:,:,gg)=sqrt(T1dFdx1.^2);
    % Un_dFdx2(:,:,gg)=sqrt(T1dFdx2.^2);
    % Un_dFdx3(:,:,gg)=sqrt(T1dFdx3.^2);
    
    % % Uncertainty in mapping function gradient due to Uai's
    % Un_dFdx1(:,:,gg)=sqrt(T2dFdx1.^2);
    % Un_dFdx2(:,:,gg)=sqrt(T2dFdx2.^2);
    % Un_dFdx3(:,:,gg)=sqrt(T2dFdx3.^2);
    
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