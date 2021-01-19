function caljob = prana_fitmodel(caljob)
%caljob = prana_fitmodel(caljob)
% wrapper function used by prana_SPIV to setup raw camera
% points identified into format expected by camera calibration
% model fitting function, and then save it back to caljob 
% for use by prana_SPIV GUI

allx1data(:,1)   = caljob.calibration_data.x_world_full{1};  % contains all x,y,z data for camera 1
allx1data(:,2)   = caljob.calibration_data.y_world_full{1};
allx1data(:,3)   = caljob.calibration_data.z_world_full{1};

allx2data(:,1)   = caljob.calibration_data.x_world_full{2};       % contains all x,y,z data for camera 2
allx2data(:,2)   = caljob.calibration_data.y_world_full{2};
allx2data(:,3)   = caljob.calibration_data.z_world_full{2};

allX1data(:,1)   = caljob.calibration_data.x_image_full{1};       % contains all X,Y data for camera 1
allX1data(:,2)   = caljob.calibration_data.y_image_full{1};

allX2data(:,1)   = caljob.calibration_data.x_image_full{2};
allX2data(:,2)   = caljob.calibration_data.y_image_full{2};

%keyboard;
% [~,in1]=sort(allx1data(:,3),1,'ascend');
% [~,in2]=sort(allx2data(:,3),1,'ascend');
% allx1data=allx1data(in1,:);
% allX1data=allX1data(in1,:);
% allx2data=allx2data(in2,:);
% allX2data=allX2data(in2,:);

% A1          = guiprops.imagemat{ind1(1)};                          % read the first image for cam 1 to get the size of the image to convert units
% A2          = guiprops.imagemat{ind2(1)};
% [rA1, cA1]  = size(A1);
% [rA2 ,cA2]  = size(A2);

rA1 = caljob.y_pixel_number;

% %JJC: I think these are unnecessary now that I treat the image coordinates
% %everywhere else in the code as pixel centered.  But is it better to just
% %fix them here and be done with it?
% allX1data(:,1)  = allX1data(:,1)-0.5;      % convert from image coords to regular coords ((0,0) at bottom left corner)
% allX1data(:,2)  = rA1-allX1data(:,2)+0.5;
% allX2data(:,1)  = allX2data(:,1)-0.5;      % convert from image coords to regular coords ((0,0) at bottom left corner)
% allX2data(:,2)  = rA1-allX2data(:,2)+0.5;

%JJC: Still need to flip the coordinate system?
%I think so since these positions are referenced to the upper left, and 
%everywhere else the image coordinates will be reference to lower left.
allX1data(:,2)  = rA1-allX1data(:,2) + 1; %plus 1 is because the first index is at 1, not 0
allX2data(:,2)  = rA1-allX2data(:,2) + 1;

modeltype   = caljob.modeltype;
optionsls   = caljob.optionsls;
optionsls   = optimset(optionsls,'Algorithm','levenberg-marquardt');  %force L-M algorithm (default is trust-region-reflective)

caljob.allx1data = allx1data;
caljob.allx2data = allx2data;
caljob.allX1data = allX1data;
caljob.allX2data = allX2data;

% fprintf('saving initial world and image coordinates....\n');
% save('calxworld.mat','allx1data');save('calyworld.mat','allx2data');save('calximage.mat','allX1data');save('calyimage.mat','allX2data');

[a_cam1, a_cam2, aXcam1, aYcam1, aXcam2, aYcam2, convergemessage] = fitcameramodels(allx1data,...
    allx2data,allX1data,allX2data,modeltype,optionsls);

%inselfcal(allx1data,allx2data,allX1data,alX2data);
caljob.aXcam1 = aXcam1;     % save the mapping coefficients
caljob.aYcam1 = aYcam1;
caljob.aXcam2 = aXcam2;
caljob.aYcam2 = aYcam2;
%[aXcam1 aYcam1 aXcam2 aYcam2]

%guiprops.scal   = 0;

guiprops.caljob.a_cam1 = a_cam1;
guiprops.caljob.a_cam2 = a_cam2;

% Do I need this section set?
%%%                                                                         %
% % zero_plane_cam1 = median(ind1);                                             %
% % zero_plane_cam2 = median(ind2);                                             %
% %                                                                             %
% % guiprops.ind1   = ind1;                                                     %
% % guiprops.ind2   = ind2;                                                     %
% % guiprops.zpc1   = zero_plane_cam1;                                          %
% % guiprops.zpc2   = zero_plane_cam2;                                          %
                                                                            %                                 %
caljob.convergemessage = convergemessage; 
