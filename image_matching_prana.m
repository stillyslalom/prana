function [deltax,deltay,NPeak,dispx,dispy,cw]=image_matching_prana(region3,region4,dx,dy)
%THIS FUNCTION CALCULATES THE DISPARITY BETWEEN TWO CLOSELY MATCHING
%IMAGES(LIKE AFTER DWO OR ITERATIVE WINDOW DEFORMATION WHEN PARTICLES
%IN IMAGE PAIR IS CLOSE TO EACH OTHER).
%INPUT:-
%REGION3=WINDOW FROM IMAGE1, WITH ZERO MEAN AND 50% GAUSSIAN SPATIAL FILTER
%OPERATION APPLIED IN PRANA.
%REGION4=WINDOW FROM IMAGE2, WITH ZERO MEAN AND 50% GAUSSIAN SPATIAL FILTER
%OPERATION APPLIED IN PRANA.
%OUTPUT:-
%DELTAX=UNCERTAINTY IN X USING DISPARITY DISTRIBUTION DISPX
%DELTAY=UNCERTAINTY IN Y USING DISPARITY DISTRIBUTION DISPY
%THIS IS BASED ON THE PRINCIPLE DISCUSSED IN "PIV UNCERTAINTY QUANTIFICATION
%BY IMAGE MATCHING", MST. 2013.BY ANDREA SCIACCHITANO,BERNHARD WIENEKE AND
%FULVIO SCARANO
%THIS CODE WRITTEN FOR PRANA IMPLEMENTATION IS WRITTEN BY SAYANTAN BHATTACHARYA.

if nargin<3
    dx = 0;
    dy = 0;
end

%% Threshold image and find peaks
sr=1; %SEARCH RADIUS FOR MATCHING PARTICLES
improduct=((region3).*(region4)); %PRODUCT OF INTENSITIES
rp3=improduct(improduct~=min(improduct(:)));%GETTING RID OF MINIMUM
if isempty(rp3)
    NPeak = 0; %all points in improduct are equal (probably 0)
else
    thresh=0.5*rms(rp3(:));%SELECTING THRESHOLD 0.5 OF RMS OF INTENSITY PRODUCT

    improduct(improduct<thresh)=0;%THRESHOLDING IMAGE

    peakmat=imregionalmax(improduct,8);%FINDING PEAKS IN 8 POINT NEIGHBOURHOOD

    cnt=peakmat(peakmat==1);
    NPeak=length(cnt);%NUMBER OF PEAKS DETECTED
end

k=1;
%% If less than 6 peaks detected reduce the threshold and iterate 3 times
while NPeak<=6 && k<=3 && ~isempty(rp3)
    
    improduct=((region3).*(region4));%PRODUCT OF INTENSITIES
    rp3=improduct(improduct~=min(improduct(:)));%GETTING RID OF MINIMUM
    if isempty(rp3)
        NPeak = 0; %all points in improduct are equal (probably 0)
    else
        thresh=(1/2^k)*rms(rp3(:));%IF NO PEAK FOUND HALF THE THRESHOLD
        improduct(improduct<thresh)=0;%THRESHOLDING IMAGE
        peakmat=imregionalmax(improduct,8);%FINDING PEAKS IN 8 POINT NEIGHBOURHOOD
        cnt=peakmat(peakmat==1);
        NPeak=length(cnt);%NUMBER OF PEAKS DETECTED
    end
    
    %fprintf('No particle found above threshold, threshold reduced to 0.5 rms');
    k=k+1;
    %return
end

% If no Peak is detected return nan values for uncertainties
if NPeak==0
    deltax=nan;
    deltay=nan;
    dispx = [];
    dispy = [];
    cw    = [];
    return;
end
    
% minimum value subtracted for 3 point gaussian fit.
region3 = region3-min(min(region3));
region4 = region4-min(min(region4));

% Find coordinates for detected peaks
[Indylist,Indxlist] =find(peakmat(:,:)==1);
%improduct2=diag(improduct(Indylist,Indxlist));% Find the intensity product
improduct2=improduct(sub2ind(size(improduct),Indylist,Indxlist));% Find the intensity product at each peak
%values at detected peak location
cw=(improduct2).^(0.5);%WEIGHTING FUNCTION

%Initialize disparity vectors
dispx=zeros(NPeak,1);
dispy=zeros(NPeak,1);


[Ny,Nx]=size(region3);% Size of window

lflag='false';% Limit flag is false when no particles are detected mear the
%edge of the window

%% FOR EACH IMAGE PRODUCT INTENSITY PEAK FIND PARTICLE CENTER LOCATIONS IN INDIVIDUAL IMAGES
for ji = 1:NPeak
    
    % Get location for this particluar peak
    Indy=Indylist(ji);
    Indx=Indxlist(ji);
    
    %IF CORNER POINTS SELECTED THEY ARE REMOVED TO AVOID TROUBLE IN
    %GAUSSIAN FITTING
    if Indy<=3 || Indy>=(Ny-3) || Indx<=3 || Indx>=(Nx-3)
        dispx(ji)=-100;
        dispy(ji)=-100;
        %cw(ji)=-100;
        %fprintf('Limits exceeded');
        lflag='true';
        %keyboard;
        continue;
    end
    
    %FORM THE WINDOW ON WHICH TO SEARCH
    sradiusy=Indy-sr:1:Indy+sr;
    sradiusx=Indx-sr:1:Indx+sr;
    
    %SELECT THE REGION FROM INPUT IMAGES WHICH CORRESPONDS TO THE WINDOW
    rA=region3(sradiusy,sradiusx);
    rB=region4(sradiusy,sradiusx);

    %Find Peaks in regions rA and rB in 4 point neighbourhood
    peakA=imregionalmax(rA,4);% 4 point neighbourhood
    peakB=imregionalmax(rB,4);% 4 point neighbourhood
    [im1y,im1x] =find(peakA==1);
    [im2y,im2x] =find(peakB==1);
    %PEAKS IN INDIVIDUAL IMAGE CORRESPONDING TO SEARCH WINDOW using maximum
%     [im1y,im1x] =find(rA==max(rA(:)));
%     [im2y,im2x] =find(rB==max(rB(:)));
    
    
    %iF MULTIPLE PEAKS FOUND SELECT THE PEAK CLOSEST TO WINDOW CENTER
    if length(im1y)>1 && length(im1x)>1
        %distance of peaks with respect to center
        d1=((im1y-(sr+1)).^2+(im1x-(sr+1)).^2).^0.5;
        in1=find(d1==min(d1));% find minimum
        in1=in1(1);
        %select minimum for window1
        im1y=im1y(in1);
        im1x=im1x(in1);
        
    end
    if length(im2y)>1 && length(im2x)>1
        %distance of peaks with respect to center
        d2=((im2y-(sr+1)).^2+(im2x-(sr+1)).^2).^0.5;
        in2=find(d2==min(d2));% find minimum
        in2=in2(1);
        %select minimum for window2
        im2y=im2y(in2);
        im2x=im2x(in2);
    end
    
    
    %GETTING PEAK COORDINATE W.R.T. IMAGE COORDINATE SYSTEM
    im1y=Indy+im1y-(sr+1);
    im1x=Indx+im1x-(sr+1);
    im2y=Indy+im2y-(sr+1);
    im2x=Indx+im2x-(sr+1);
    
    %DOING SUBPIXEL FIT FOR EACH PEAK EACH DIRECTION
    % for ith peak in image1, subpixel fit in y direction
    lCm11 = log(region3(im1y-1,im1x));
    lC001 = log(region3(im1y,im1x));
    lCp11 = log(region3(im1y+1,im1x));
    if (2*(lCm11+lCp11-2*lC001)) == 0 % if three point Gaussian fit fail, set the particle diameter to a deafault guess value
        shift_erry1 = 0;
    else
        shift_erry1 = (lCm11-lCp11)/(2*(lCm11+lCp11-2*lC001));
    end
    clear lCm11 lC001 lCp11;
    % for ith peak in image1, subpixel fit in x direction
    lCm11 = log(region3(im1y,im1x-1));
    lC001 = log(region3(im1y,im1x));
    lCp11 = log(region3(im1y,im1x+1));
    if (2*(lCm11+lCp11-2*lC001)) == 0 % if three point Gaussian fit fail, set the particle diameter to a deafault guess value
        shift_errx1 = 0;
    else
        shift_errx1 = (lCm11-lCp11)/(2*(lCm11+lCp11-2*lC001));
    end
    clear lCm11 lC001 lCp11;
    % for ith peak in image2, subpixel fit in y direction
    lCm11 = log(region4(im2y-1,im2x));
    lC001 = log(region4(im2y,im2x));
    lCp11 = log(region4(im2y+1,im2x));
    if (2*(lCm11+lCp11-2*lC001)) == 0 % if three point Gaussian fit fail, set the particle diameter to a deafault guess value
        shift_erry2 = 0;
    else
        shift_erry2 = (lCm11-lCp11)/(2*(lCm11+lCp11-2*lC001));
    end
    clear lCm11 lC001 lCp11;
    % for ith peak in image2, subpixel fit in x direction
    lCm11 = log(region4(im2y,im2x-1));
    lC001 = log(region4(im2y,im2x));
    lCp11 = log(region4(im2y,im2x+1));
    if (2*(lCm11+lCp11-2*lC001)) == 0 % if three point Gaussian fit fail, set the particle diameter to a deafault guess value
        shift_errx2 = 0;
    else
        shift_errx2 = (lCm11-lCp11)/(2*(lCm11+lCp11-2*lC001));
    end
    clear lCm11 lC001 lCp11;
    
    
    %PARTICLE CENTER LOCATIONS UPTO SUBPIXEL APPROXIMATION
    Indy1=im1y+shift_erry1;
    Indx1=im1x+shift_errx1;
    Indy2=im2y+shift_erry2;
    Indx2=im2x+shift_errx2;
    
    % If subpixel fit is nan or greater than 1 then subpixel fit failed,
    % that point neglected
    if abs(shift_errx1)>1 || abs(shift_errx2)>1 || abs(shift_erry1)>1 || abs(shift_erry2)>1
        dispx(ji)=-200;
        dispy(ji)=-200;
    elseif isnan(shift_errx1) || isnan(shift_errx2) || isnan(shift_erry1) || isnan(shift_erry2)
        dispx(ji)=-300;
        dispy(ji)=-300;
    else
        %DISPARITY BETWEEN PARTICLE PAIR FOR EACH PAIR OF MATCHED PARTICLE
        dispx(ji)=Indx2-Indx1;
        dispy(ji)=Indy2-Indy1;
    end
    
end

%% DISCARDING POINTS WHICH WERE CLOSE TO THE BOUNDARY
if strcmp(lflag,'true')
    %keyboard;
    ind=find(dispx(:)==-100);
    dispx(ind)=[];
    dispy(ind)=[];
    cw(ind)=[];
    %NPeak
    NPeak=NPeak-length(ind);
    %length(ind)
    %NPeak
end
%% DISCARDING POINTS WHERE GAUSSIAN FIT FAILED
ind2=find(dispx(:)==-200);
dispx(ind2)=[];
dispy(ind2)=[];
cw(ind2)=[];
%NPeak
NPeak=NPeak-length(ind2);

ind3=find(dispy(:)==-300);
dispx(ind3)=[];
dispy(ind3)=[];
cw(ind3)=[];
%NPeak
NPeak=NPeak-length(ind3);

%% If after this elimination no Peak is detected return nan values for uncertainties
if NPeak==0
    deltax=nan;
    deltay=nan;
    return;
end

%% Remove any subpixel shift left in the images
dispx = dispx - dx;
dispy = dispy - dy;

%% CALCULATING WEIGHTED MEAN AND VARIANCE OF DISPARITY DISTRIBUTION FOR THIS
%WINDOW PAIRS AND SAVING THE UNCERTAINTY AS OUTPUT.
mewx=(1/sum(cw))*sum(cw.*dispx);% Weiighted mean
sigx=sqrt((1/sum(cw)).*sum(cw.*((dispx-mewx).^2)));% Weighted std deviation
deltax=sqrt(mewx^2 + (sigx/sqrt(NPeak))^2);% Uncertainty in x

mewy=(1/sum(cw))*sum(cw.*dispy);% Weiighted mean
sigy=sqrt((1/sum(cw)).*sum(cw.*((dispy-mewy).^2)));% Weighted std deviation
deltay=sqrt(mewy^2 + (sigy/sqrt(NPeak))^2);% Uncertainty in y

end