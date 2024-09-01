%% Import options
cd '\\131.155.50.88\Larry Fitzpatrick\20231216_MCF7_4-colour_EGF-PDL1-Combo-Tf\Analysis\Tf_Analysis'   % Set directory
N=14;    % Set number of folders with images (roisets)
filenameSpots='Tf_%d.csv';    % Set name of trajectory file (change the for %d). They should be numbered from 1 to N
% filenameStat='Track statistics %d.2.csv';
% namePCA = 'PCA_A549.csv';

% Analysis options
datasetCode = 3;    % Add a numeric code for each sample types (e.g. each cell line), to be added in the first column of the single-cell matrix
nPar = 14;          % Number of total analysis paramets (do not change this number if the code below is not modified)

% Variables
pixelsize = 117;    % Size of pixel in nm
expT = 0.05;        % Exposure time in s
cx = 5;             % Column of x positions in file spots
cy = 6;             % Column of y positions in file spots
ci = 3;             % Column of track index in file spots

% Output data saving names
nameTrackInfoRed = 'ROI_a_%d_Track_Info_Red.csv';         % Tracks information of a single ROI
nameTrackInfoAllRed = 'Track_a_Information_All_Red.csv';  % Combined tack information of all ROIs
nameCellInfoMeanRed = 'Cell_a_Information_All_Red.csv';   % Average single cell infromation (includes track information, plus density and step count)
nameCellInfoMedianRed = 'Cell_a_Information_Median_Red.csv';  
nameCellInfoSTDRed = 'Cell_a_Information_STD_Red.csv';  

nameTrackInfoGreen = 'ROI_a_%d_Track_Info_Green.csv';         % Tracks information of a single ROI
nameTrackInfoAllGreen = 'Track_a_Information_All_Green.csv';  % Combined tack information of all ROIs
nameCellInfoMeanGreen = 'Cell_a_Information_All_Green.csv';   % Average single cell infromation (includes track information, plus density and step count)
nameCellInfoMedianGreen = 'Cell_a_Information_Median_Green.csv';  
nameCellInfoSTDGreen = 'Cell_a_Information_STD_Green.csv';  

%% %%%%%%%%%%%%%%%%%%% ROI SLECTION %%%%%%%%%%%%%%%%%%
zipname = 'RoiSet_%d';
NameROIgreen='ROIspotsgreen_%d.csv';
NameROIred='ROIspotsred_%d.csv';
% NameStatgreen='ROIstatsgreen_%d.csv';
% NameStatred='ROIstatsred_%d.csv';
CountROI=0;
area=[];

for i = 1:N
%% Data Import
fileSpots = sprintf(filenameSpots,i);
data = readmatrix(fileSpots);
if numel(unique(data(:,1)))==1  % If there is only one channel
           dataRed=data;
           
           n=length(dataRed);
           A=zeros(n,10);
           A(:,1)=dataRed(:,1);
           A(:,2)=dataRed(:,2);
           A(:,3)=dataRed(:,3);
           A(:,4)=dataRed(:,4);
           A(:,5)=dataRed(:,5);
           A(:,6)=dataRed(:,6);
           A(:,7)=dataRed(:,7);
           A(:,8)=dataRed(:,8);
           A(:,9)=dataRed(:,9);
           A(:,10)=dataRed(:,10);
         

        else       % If there are more than one channel
            
           disp('More than one channel!');
           dataGreen=data(data(:,1)==0,:); % Green
           dataRed=data(data(:,1)==1,:);; %Red
           
            n=length(dataGreen);
            nref=length(dataRed);
            A=zeros(n,10);
            B=zeros(nref,10);
            A(:,1)=dataGreen(:,1);
            A(:,2)=dataGreen(:,2);
            A(:,3)=dataGreen(:,3);
            A(:,4)=dataGreen(:,4);
            A(:,5)=dataGreen(:,5);
            A(:,6)=dataGreen(:,6);
            A(:,7)=dataGreen(:,7);
            A(:,8)=dataGreen(:,8);
            A(:,9)=dataGreen(:,9);
            A(:,10)=dataGreen(:,10);
            B(:,1)=dataRed(:,1);
            B(:,2)=dataRed(:,2);
            B(:,3)=dataRed(:,3);
            B(:,4)=dataRed(:,4);
            B(:,5)=dataRed(:,5);
            B(:,6)=dataRed(:,6);
            B(:,7)=dataRed(:,7);
            B(:,8)=dataRed(:,8);
            B(:,9)=dataRed(:,9);
            B(:,10)=dataRed(:,10);
end
            
%  trnumgreen=A(:,3);
%  trnumgreen(:,2)=0;
%  trposgreen=zeros(length(trnumgreen),2);
%  for n=1:length(trnumgreen)
%     trnumgreen(n,2)=find(data(:,3)==trnumgreen(n,1),1,'first');
%     trposgreen(n,1)=data(trnumgreen(n,2),5);
%     trposgreen(n,2)=data(trnumgreen(n,2),6);
%  end
%  trposgreen=trposgreen/117;

trnumred=A(:,3);
trnumred(:,2)=0;
trposred=zeros(length(trnumred),2);
for n=1:length(trnumred)
   trnumred(n,2)=find(data(:,3)==trnumred(n,1), 1);
   trposred(n,1)=data(trnumred(n,2),5);
   trposred(n,2)=data(trnumred(n,2),6);
end
trposred=trposred/117;
% fileStatgreen = sprintf(filenameStatgreen,i);
% fileStatred = sprintf(filenameStatred,i);
% dataStatgreen = readmatrix(fileStatgreen); 
% dataStatred = readmatrix(fileStatred); 

%% Extract Positions from ROIs
    zipfile=sprintf(zipname,i);
    zipfilefolder=strcat(cd,'/',zipfile);
    addpath(zipfilefolder);
    
    % Obtain Positions From ROIs and save them in CSV files
    Fext = dir([zipfilefolder '/' '*ext.mat']);
    Fext = struct2table(Fext);
    Fint = dir([zipfilefolder '/' '*int.mat']);
    Fint = struct2table(Fint);
    for k = 1:size(Fext,1)
        CountROI=CountROI+1;
%         check_roigreen=zeros(size(A,1),1); 
        check_roired=zeros(size(A,1),1);
        load(char(Fext(k,1).name));
        % Crop inner roi
%          r = find(strcmp(Fint.name,sprintf('roi_%d_int.mat',k)));
%          if r == true
%             load(char(Fint(r,1).name));
%             rmi = rmi*-1;
%             rme = rme + rmi;
%          end
        % Select data in matrix
%          for n=1:size(trposgreen,1)
%              check_roigreen(n,1)=rme(fix(trposgreen(n,2))+1,fix(trposgreen(n,1))+1);
%          end
        for n=1:size(trposred,1)
            check_roired(n,1)=rme(fix(trposred(n,2))+1,fix(trposred(n,1))+1);
        end
        area(CountROI,1)=sum(rme(:)==1)*((pixelsize/1000)^2);   % Save area of each ROI in um2
        % Obtain ROI Localizations and save them
%         dataROIgreen=A(any(check_roigreen(:,1)==1,2),:);
        dataROIred=A(any(check_roired(:,1)==1,2),:);
%         csvwrite(sprintf(NameROIgreen,CountROI),dataROIgreen);
        csvwrite(sprintf(NameROIred,CountROI),dataROIred);
        
%         % Save trajectory statistics
%          index = unique(dataROIgreen(:,ci));
%          statROIgreen = [];
%          for m = 1:size(index,1)
%             statROIgreen(m,:) = dataStatgreen(find(dataStatgreen(:,2)==index(m)),:); 
%          end
%          csvwrite(sprintf(NameStatgreen,CountROI),statROI);
%          
%          index = unique(dataROIred(:,ci));
%          statROIred = [];
%          for m = 1:size(index,1)
%             statROIred(m,:) = dataStatred(find(dataStatred(:,2)==index(m)),:); 
%          end
%          csvwrite(sprintf(NameStat,CountROI),statROI);
    end
end

%% Data analysis and gathering (RED)
% Initialize variables and arrays
cellInfoMean_Red = zeros(CountROI, nPar+3);     % Make matrix to store average cell information (with 3 extra column for ID, density and steps)
cellInfoMean_Red(:,1) = datasetCode;            % Add dataset code to first column
trackInfoAll_Red = [];                          % Initialize matrix to save track data from all images
cellInfoMedian_Red = zeros(CountROI, nPar+1);     % Make matrix to store average cell information (with 3 extra column for ID, density and steps)
cellInfoMedian_Red(:,1) = datasetCode;            % Add dataset code to first column
cellInfoSTD_Red = zeros(CountROI, nPar+1);     % Make matrix to store average cell information (with 3 extra column for ID, density and steps)
cellInfoSTD_Red(:,1) = datasetCode;            % Add dataset code to first column


for i = 1:CountROI
    % Opend files
    fileROIred = sprintf(NameROIred, i);
    dataROIred = readmatrix(fileROIred);
    if ~isempty(dataROIred)
        dataROIred(:, cx:cy) = dataROIred(:, cx:cy) / 1000;  % Transform nm to um of X and Y positions
        % Initialize file where data would be sotred
        trackID = unique(dataROIred(:,ci));
        trackInfo = zeros(size(trackID,1), nPar+1);
        trackInfo(:,1) = i;          % Add ROI ID to first column

        % Loop around each track
        for m = 1:size(trackID,1)
            % Get data from track
            trackData = dataROIred(dataROIred(:,ci)==trackID(m,1),:);

            % Track displacement, um (column 2)
            trackInfo(m,2) = pdist2([trackData(1,cx) , trackData(1,cy)] , [trackData(end,cx) , trackData(end,cy)]);

            % Track duration, s (colum 3)
            trackInfo(m,3) = (trackData(end,2) - trackData(1,2)+1) * expT; 

            % Loop the steps per track
            stepInfo = zeros(size(trackData,1)-1 , 5);    % Temporary file to store data
            setpInfo(1,1) = NaN;                        % First angles is not counted

            for n = 1:size(trackData,1)-1   
                % Measure angles (there is one less angles than steps, so first one is skipped)
                if n > 1
                   vec1 = [trackData(n-1,cx)-trackData(n,cx) , trackData(n-1,cy)-trackData(n,cy) , 0];
                   vec2 = [trackData(n+1,cx)-trackData(n,cx) , trackData(n+1,cy)-trackData(n,cy) , 0];
                   stepInfo(n,1) = atan2(norm(cross(vec1,vec2)),dot(vec1,vec2));
                end

                stepInfo(n,2) = pdist2([trackData(n+1,cx) , trackData(n+1,cy)] , [trackData(n,cx) , trackData(n,cy)]);  % For Total Distance Travelled
                stepInfo(n,3) = pdist2([trackData(n+1,cx) , trackData(n+1,cy)] , [trackData(1,cx) , trackData(1,cy)]);  % For Max Distance Travelled
                stepInfo(n,4) = stepInfo(n,2) / ((trackData(n+1,2) - trackData(n,2)) * expT);                           % For Track Speeds
                stepInfo(n,5) = (trackData(n+1,cx) - trackData(n,cx))^2 + (trackData(n+1,cy) - trackData(n,cy))^2;      % For Mean Square Displacement and Diffusion Coefficient       
            end
            % Mean Directionality Change Rate, rad (column 4)
            trackInfo(m,4) = nanmean(stepInfo(:,1),'all');
            % Total Distance Traveled, um (column 5)
            trackInfo(m,5) = sum(stepInfo(:,2));
            % Max Distance Travelled, um (column 6)
            trackInfo(m,6) = max(stepInfo(:,3));
            % Track Speeds, um/s
            trackInfo(m,7) = mean(stepInfo(:,4));       % Mean track speed (column 7)
            trackInfo(m,8) = min(stepInfo(:,3));        % Minimum track speed (column 8)
            trackInfo(m,9) = max(stepInfo(:,3));        % Maximum track speed (column 9)
            trackInfo(m,10) = median(stepInfo(:,3));    % Median track speed (column 10)
            % Average Mean Square displacement and average instant Diffusion Coefficient
            trackInfo(m,11) = mean(stepInfo(:,5));      % MSD, um2(column 11)
            trackInfo(m,12) = trackInfo(m,11)/4*expT;   % Diffusion Coeeficient, um2/s (column 12)
            % Confinement Ratio (column 13)
            trackInfo(m,13) = trackInfo(m,2) / trackInfo(m,5);
            % Mean Straight Line Speed, um/s (column 14)
            trackInfo(m,14) = trackInfo(m,3) / trackInfo(m,2);
            % Linearity of forward Progression (column 15)
            trackInfo(m,15) = trackInfo(m,14) / trackInfo(m,7);
        end
    % Save Track Info ROI
    saveTrackInfo = sprintf(nameTrackInfoRed , i);
    csvwrite(saveTrackInfo , trackInfo);

    % Append Track infomration
    trackInfoAll_Red = [trackInfoAll_Red ; trackInfo];

    else
        % If dataROIgreen is empty, save an empty matrix (empty files are not appended in to the global matrix)
        trackInfo = NaN;
        saveTrackInfo = sprintf(nameTrackInfoRed , i);
        csvwrite(saveTrackInfo , trackInfo);
    end
    
    %% Average values per Cell
    if ~isempty(dataROIred)
        for n= 2:nPar+1
            cellInfoMean_Red(i,n) = mean(trackInfo(:,n));    
        end
        % Add Density (column 16) and mean steps per cell (column 17)
        % Calculate Density     
        cellInfoMean_Red(i,16) = numel(unique(dataROIred(:, ci))) / area(i, 1);
%         % Calculate number of steps per cell
%         cellInfoMean_Red(i,17) = numel(dataROIred(:, ci)) - numel(unique(dataROIred(:, ci)));

    else
        cellInfoMean_Red(i,:) = [];  % In case dataset is empty, remove line from matrix
    end
    %% Median values per cell 
    if ~isempty(dataROIred)
        for n= 2:nPar+1
            cellInfoMedian_Red(i,n) = median(trackInfo(:,n)); 
        end       
 
    else
        cellInfoMedian_Red(i,:) = [];  % In case dataset is empty, remove line from matrix
    end
    
    %% STD values per cell
    if ~isempty(dataROIred)
        for n= 2:nPar+1
            cellInfoSTD_Red(i,n) = std(trackInfo(:,n)); 
        end       
 
    else
        cellInfoSTD_Red(i,:) = [];  % In case dataset is empty, remove line from matrix
    end
end

% Save Track Information All
csvwrite(nameTrackInfoAllRed , trackInfoAll_Red);
csvwrite(nameCellInfoMeanRed , cellInfoMean_Red);
csvwrite(nameCellInfoMedianRed, cellInfoMedian_Red);
csvwrite(nameCellInfoSTDRed, cellInfoSTD_Red);


% %% Data analysis and gathering (GREEN)
% % Initialize variables and arrays
% cellInfoMean_Green = zeros(CountROI, nPar+3);     % Make matrix to store average cell information (with 3 extra column for ID, density and steps)
% cellInfoMean_Green(:,1) = datasetCode;            % Add dataset code to first column
% trackInfoAll_Green = [];                          % Initialize matrix to save track data from all images
% cellInfoMedian_Green = zeros(CountROI, nPar+1);     % Make matrix to store average cell information (with 3 extra column for ID, density and steps)
% cellInfoMedian_Green(:,1) = datasetCode;            % Add dataset code to first column
% cellInfoSTD_Green = zeros(CountROI, nPar+1);     % Make matrix to store average cell information (with 3 extra column for ID, density and steps)
% cellInfoSTD_Green(:,1) = datasetCode;            % Add dataset code to first column
% 
% 
% for i = 1:CountROI
%     % Opend files
%     fileROIgreen = sprintf(NameROIgreen, i);
%     dataROIgreen = readmatrix(fileROIgreen);
%     if ~isempty(dataROIgreen)
%         dataROIgreen(:, cx:cy) = dataROIgreen(:, cx:cy) / 1000;  % Transform nm to um of X and Y positions
%         % Initialize file where data would be sotred
%         trackID = unique(dataROIgreen(:,ci));
%         trackInfo = zeros(size(trackID,1), nPar+1);
%         trackInfo(:,1) = i;          % Add ROI ID to first column
% 
%         % Loop around each track
%         for m = 1:size(trackID,1)
%             % Get data from track
%             trackData = dataROIgreen(dataROIgreen(:,ci)==trackID(m,1),:);
% 
%             % Track displacement, um (column 2)
%             trackInfo(m,2) = pdist2([trackData(1,cx) , trackData(1,cy)] , [trackData(end,cx) , trackData(end,cy)]);
% 
%             % Track duration, s (colum 3)
%             trackInfo(m,3) = (trackData(end,2) - trackData(1,2)+1) * expT; 
% 
%             % Loop the steps per track
%             stepInfo = zeros(size(trackData,1)-1 , 5);    % Temporary file to store data
%             setpInfo(1,1) = NaN;                        % First angles is not counted
% 
%             for n = 1:size(trackData,1)-1   
%                 % Measure angles (there is one less angles than steps, so first one is skipped)
%                 if n > 1
%                    vec1 = [trackData(n-1,cx)-trackData(n,cx) , trackData(n-1,cy)-trackData(n,cy) , 0];
%                    vec2 = [trackData(n+1,cx)-trackData(n,cx) , trackData(n+1,cy)-trackData(n,cy) , 0];
%                    stepInfo(n,1) = atan2(norm(cross(vec1,vec2)),dot(vec1,vec2));
%                 end
% 
%                 stepInfo(n,2) = pdist2([trackData(n+1,cx) , trackData(n+1,cy)] , [trackData(n,cx) , trackData(n,cy)]);  % For Total Distance Travelled
%                 stepInfo(n,3) = pdist2([trackData(n+1,cx) , trackData(n+1,cy)] , [trackData(1,cx) , trackData(1,cy)]);  % For Max Distance Travelled
%                 stepInfo(n,4) = stepInfo(n,2) / ((trackData(n+1,2) - trackData(n,2)) * expT);                           % For Track Speeds
%                 stepInfo(n,5) = (trackData(n+1,cx) - trackData(n,cx))^2 + (trackData(n+1,cy) - trackData(n,cy))^2;      % For Mean Square Displacement and Diffusion Coefficient       
%             end
%             % Mean Directionality Change Rate, rad (column 4)
%             trackInfo(m,4) = nanmean(stepInfo(:,1),'all');
%             % Total Distance Traveled, um (column 5)
%             trackInfo(m,5) = sum(stepInfo(:,2));
%             % Max Distance Travelled, um (column 6)
%             trackInfo(m,6) = max(stepInfo(:,3));
%             % Track Speeds, um/s
%             trackInfo(m,7) = mean(stepInfo(:,4));       % Mean track speed (column 7)
%             trackInfo(m,8) = min(stepInfo(:,3));        % Minimum track speed (column 8)
%             trackInfo(m,9) = max(stepInfo(:,3));        % Maximum track speed (column 9)
%             trackInfo(m,10) = median(stepInfo(:,3));    % Median track speed (column 10)
%             % Average Mean Square displacement and average instant Diffusion Coefficient
%             trackInfo(m,11) = mean(stepInfo(:,5));      % MSD, um2(column 11)
%             trackInfo(m,12) = trackInfo(m,11)/4*expT;   % Diffusion Coeeficient, um2/s (column 12)
%             % Confinement Ratio (column 13)
%             trackInfo(m,13) = trackInfo(m,2) / trackInfo(m,5);
%             % Mean Straight Line Speed, um/s (column 14)
%             trackInfo(m,14) = trackInfo(m,3) / trackInfo(m,2);
%             % Linearity of forward Progression (column 15)
%             trackInfo(m,15) = trackInfo(m,14) / trackInfo(m,7);
%         end
%     % Save Track Info ROI
%     saveTrackInfo = sprintf(nameTrackInfoGreen , i);
%     csvwrite(saveTrackInfo , trackInfo);
% 
%     % Append Track infomration
%     trackInfoAll_Green = [trackInfoAll_Green ; trackInfo];
% 
%     else
%         % If dataROIgreen is empty, save an empty matrix(empty files are not appended in to the global matrix)
%         trackInfo = NaN;
%         saveTrackInfo = sprintf(nameTrackInfoGreen , i);
%         csvwrite(saveTrackInfo , trackInfo);
%     end
%     
%     %% Average values per Cell
%     if ~isempty(dataROIgreen)
%         for n= 2:nPar+1
%             cellInfoMean_Green(i,n) = mean(trackInfo(:,n)); 
%         end
%         % Add Density (column 16) and mean steps per cell (column 17)
%         % Calculate Density     
%         cellInfoMean_Green(i,16) = numel(unique(dataROIgreen(:, ci))) / area(i, 1);
% %         % Calculate number of steps per cell
% %         cellInfoMean_Green(i,17) = numel(dataROIgreen(:, ci)) - numel(unique(dataROIgreen(:, ci)));
%         
%  
%     else
%         cellInfoMean_Green(i,:) = [];  % In case dataset is empty, remove line from matrix
%     end
%     
%     %% Median values per cell 
%     if ~isempty(dataROIgreen)
%         for n= 2:nPar+1
%             cellInfoMedian_Green(i,n) = median(trackInfo(:,n)); 
%         end       
%  
%     else
%         cellInfoMedian_Green(i,:) = [];  % In case dataset is empty, remove line from matrix
%     end
%     
%     %% STD values per cell
%     if ~isempty(dataROIgreen)
%         for n= 2:nPar+1
%             cellInfoSTD_Green(i,n) = std(trackInfo(:,n)); 
%         end       
%  
%     else
%         cellInfoSTD_Green(i,:) = [];  % In case dataset is empty, remove line from matrix
%     end
% end
% 
% % Save Track Information All
% csvwrite(nameTrackInfoAllGreen , trackInfoAll_Green);
% csvwrite(nameCellInfoMeanGreen , cellInfoMean_Green);
% csvwrite(nameCellInfoMedianGreen, cellInfoMedian_Green);
% csvwrite(nameCellInfoSTDGreen, cellInfoSTD_Green);
