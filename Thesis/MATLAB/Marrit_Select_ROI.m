close all
%% Importing options
cd '\\131.155.50.88\Larry Fitzpatrick\20231216_MCF7_4-colour_EGF-PDL1-Combo-Tf\B10-B11_EGFR_PDL1_1nM\pos_5'    % Set directory
BF = 'MCF7_4-colour_EGF-PDL1-Combo-Tf_B10_EGFR_1nM_posXY5_channels_t0_posZ0_colour0.tif';    % Set name of brigthfield image
BF2 = 'MCF7_4-colour_EGF-PDL1-Combo-Tf_B10_EGFR_1nM_posXY5_channels_t0_posZ0_colour1_thumbnail.png';
filename = 'EGFR_PDL1_05.csv';    % Set name of localizations / trajectories 
cx = 5;  % Column number of x position in file
cy = 6;  % Column number of y position in file
pixelsize = 117;    % Pixel size in nm

BaseFolder = '\\131.155.50.88\Larry Fitzpatrick\20231216_MCF7_4-colour_EGF-PDL1-Combo-Tf\B10-B11_EGFR_PDL1_1nM\pos_5';

%% Import Localizations and BF
%data = readmatrix(filename); 
%data(:,cx)=data(:,cx);        % Transform nm into pixel
%data(:,cy)=data(:,cy);
[X,cmap] = imread(BF);
[X2, cmap2] = imread(BF2);
TR=csvread(filename,1,0);
trnum=unique(TR(:,3));
trnum(:,2)=0;
trpos=zeros(length(trnum),2);
for i=1:length(trnum)
   trnum(i,2)=find(TR(:,3)==trnum(i,1), 1, 'first');
   trpos(i,1)=TR(trnum(i,2),5);
   trpos(i,2)=TR(trnum(i,2),6);
end
trpos=trpos/117;

%% Plot Image
figure
imshow(X,cmap); hold on
scatter(trpos(:,1),trpos(:,2),10,'r','filled','o');

%figure
%imshow(X,cmap); hold on
%scatter(data(:,cx),data(:,cy),10,'r','filled','o');

%figure
%imshow(X2,cmap2); hold on
%scatter(data(:,cx),data(:,cy),10,'r','filled','o');
