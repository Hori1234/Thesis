
%% Prepare folder (run this section changing the ROISet number for each image)
mkdir ROISet_04;
BaseFolder=strcat(BaseFolder,'/ROISet_04');
addpath(BaseFolder)

%% Crop and save (copy and run one of the options to make a ROI)
crop=imellipse
crop=imfreehand
crop=impoly
crop=imrect

%% Save Area's as appropriate
rme=createMask(crop); % External
%rmi=createMask(crop); % Internal

%% Save cropping (run this section to save the roi, changing the roi number, starting from 1 for each image)

save(strcat(BaseFolder,'/roi_04_ext'),'rme'); % External

% If there is an internal roi, run this part with the same number as the external roi
%save(strcat(BaseFolder,'/roi_14_int'),'rmi');