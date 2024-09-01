%% Import options

cd '\\131.155.50.88\Larry Fitzpatrick\20231216_MCF7_4-colour_EGF-PDL1-Combo-Tf\Analysis\EGFR_Analysis'   % Set directory: special analysis folder:just numbered ROI sets and corresp. numbered track files (=csv from ONI)
N=14;    % Set number of folders with images (roisets)
filenameSpots='EGFR_%d.csv';    % Set name of trajectory file(=csv from ONI) (change what is before the %d). Numbered from 1 to N
% filenameStat='Track statistics %d.2.csv';
% namePCA = 'PCA_A549.csv';%not used

% Analysis options
minStep = 3;    % Minimum number of steps or frames to do be considered right 
cellIdx = 4; % 1 for MDA468, 2 for MDA231 and 3 for MCF7 4 for A431 5 for RPTEC (Larry) 6 for U251MG (Larry)

pixelsize = 117;        % Size of pixel in nm careful cause also hardcoded later on 
FOV=[428 500];          % FOV size in pixels (you can check brigthfield images imageJ)
expT = 0.03;    % Exposure time in s
cx = 5;     % Column of x positions in file spots
cy = 6;     % Column of y positions in file spots
ci = 3;     % Column of track index in file spots

%% %%%%%%%%%%%%%%%%%%% ROI SLECTION %%%%%%%%%%%%%%%%%%
zipname = 'ROISet_%d';%name your folders with the ROIs for each FOV this
NameROIgreen='ROIspotsint_Tf_col_0_Full_%d.csv';
NameROIred='ROIspotsint_PDL1_col_1_Full_%d.csv';
% NameStatgreen='ROIstatsgreen_%d.csv';
% NameStatred='ROIstatsred_%d.csv';
CountROI=0;
area=[];


for i = 1:N%repeat for each RoiSet
%% Data Import
fileSpots = sprintf(filenameSpots,i);%sprintf=string print formatted:replaces placeholder %d (=integers) in filenameSpots with value of i
data = readmatrix(fileSpots);%read the ONI csv file and put into data
if numel(unique(data(:,1)))==1  % If there is only one channel
    %numel returns nr of elements in array, unique returns all unique
    %values so if the nr of unique values in channelID is one
           dataRed=data;
           A=dataRed;
           A(:,11)=[];%remove 11th column 
          else       % If there are 2 channels
           disp('More than one channel!');
           dataGreen=data(data(:,1)==0,:);%put all data with channelID zero into dataGREEN
           dataRed=data(data(:,1)==1,:); %same for red
           A=dataGreen;
           A(:,11)=[];%remove 11th column
           B=dataRed;
           B(:,11)=[];
end
  %store green X and Y coordinates and create A2 with all columns for unique red tracks only          
 trnumgreen=A(:,3);%store column 3 (trackIDdata)
 trnumgreenU=unique(trnumgreen); % store only unique trackIDs
 A2=zeros(length(trnumgreenU),11);%new matrix filled with zeros to store all data from unique tracks only
 trposgreenU=zeros(length(trnumgreenU),2);
 trnumgreenU(:,2)=0; %add second column filled with zeros
 
 for n=1:length(trnumgreenU)%so for all unique green tracks
    trnumgreenU(n,2)=find(dataGreen(:,3)==trnumgreenU(n,1),1,'first');%find the INDEX of first time trackID column data=same as nth in trnumgreen
    %store for now in 2nd colum 'first'element of index i.e the row nr. 
    trposgreenU(n,1)=dataGreen(trnumgreenU(n,2),cx);%get Xcoordinate from 5th column using rowNr stored above
    A2(n,:)=dataGreen(trnumgreenU(n,2),:);%store all in A2 for using rowNR stored above
    trposgreenU(n,2)=dataGreen(trnumgreenU(n,2),cy);% replace rowNR with Y coordinate
 end
 trposgreenU=trposgreenU/pixelsize;%divide by pixelsize to convert to nm
A2(:,11)=[];%also remove 11th column in table/matrix containing only unique tracks
 
%store red X and Y coordinates and create B2 with all columns for unique red tracks only   
 trnumred=B(:,3);%exactly same for red
 trnumredU=unique(trnumred); 
 B2=zeros(length(trnumredU),11);
 trposredU=zeros(length(trnumredU),2);
 trnumredU(:,2)=0;

 for n=1:length(trnumredU)
    trnumredU(n,2)=find(dataRed(:,3)==trnumredU(n,1),1,'first');
    trposredU(n,1)=dataRed(trnumredU(n,2),cx); 
    B2(n,:)=dataRed(trnumredU(n,2),:); 
    trposredU(n,2)=dataRed(trnumredU(n,2),cy);
 end
B2(:,11)=[];
trposredU=trposredU/pixelsize;
% fileStatgreen = sprintf(filenameStatgreen,i);
% fileStatred = sprintf(filenameStatred,i);
% dataStatgreen = readmatrix(fileStatgreen); 
% dataStatred = readmatrix(fileStatred); 
%% Extract Positions from ROIs
    zipfile=sprintf(zipname,i);%replace %d in zipname with real value (i)
    zipfilefolder=strcat(cd,'/',zipfile);%STRingconCATenate so use to open right zipfilefolder containing ROI
    addpath(zipfilefolder);
    
    % Obtain Positions From ROIs and save them in CSV files
    Fext = dir([zipfilefolder '/' '*ext.mat']);%dir lists files in directory lists files (because of slash added in dir looking inside) zipfilefolder that end in int.mat
    %gives structure matrix with name, folder and more info
    Fext = struct2table(Fext); %convert structure matrix into table
    Fint = dir([zipfilefolder '/' '*int.mat']);
    Fint = struct2table(Fint);
    for k = 1:size(Fext,1)% fancy way to determine nr Rows in Fext because is table so for each ext.mat file in the folder
        CountROI=CountROI+1;%add 1 to CountROi for each iteration loop
        check_roigreenU=zeros(size(trnumgreenU,1),1); 
        check_roiredU=zeros(size(trnumredU,1),1);
        load(char(Fext(k,1).name)); %from the filename column ('name') of Fext extract kth row as a character array and load that:
        %gives rme: 500x428 (=pixels in image) logical: each pixel 1 if in ROI 0 if not
        % Crop inner roi
         r = find(strcmp(Fint.name,sprintf('roi_%d_int.mat',k)));%strcmp compares strings so save index where kth int.mat is 
         if r == true %so if r=1 meaning if there is an inner ROI
         load(char(Fint(r,1).name));%load int roi
         rmi = rmi*-1; %multiply the logical rmi with -1
         rme = rme + rmi;%add the now neg. rmi to rme so wherever rme was 1 but rmi was also 1 (and now -1) rme becomes 0 
        end
        % Select data in matrix
         for n=1:size(trnumgreenU,1)%for all unique trackIDs 
            check_roigreenU(n,1)=rme(fix(trposgreenU(n,2))+1,fix(trposgreenU(n,1))+1);%find XY cooridnates this trackID in rme and 
            %only when result is 1 (so inside outer ROI and not in inner ROI) save coordinates in check_roigreen
         end
        for n=1:size(trnumredU,1) %same for red
            %weird:had to switch X and Y this step: reliable?
           check_roiredU(n,1)=rme(fix(trposredU(n,2))+1,fix(trposredU(n,1))+1);
        end
        area(CountROI,1)=sum(rme(:)==1)*((pixelsize/1000)^2);   % Save area of each ROI in um2
        % Obtain ROI Localizations and save them
        dataROIgreenU=A2(any(check_roigreenU(:,1)==1,2),:);%use any to find wherever check_roigreen is not 0:
        %only then copy all columns from A2 for this trackID and save
        dataROIredU=B2(any(check_roiredU(:,1)==1,2),:); 
        csvwrite(sprintf(NameROIgreen,CountROI),dataROIgreenU);%now save file as cvs
        csvwrite(sprintf(NameROIred,CountROI),dataROIredU);
        
        % Save trajectory statistics
%         index = unique(dataROI(:,ci));
%         statROI = [];
%         for m = 1:size(index,1)
%            statROI(m,:) = dataStat(find(dataStat(:,2)==index(m)),:); 
%         end
%         csvwrite(sprintf(NameStat,CountROI),statROI);
    end
end

%% %%%%%%%%%%%%%%%%%%%%% Data Gathering %%%%%%%%%%%%%%%%%%%%%%%%%%%
% dataPCA = []; 
% dataPCA(:,1) = cellIdx;
% density = zeros(CountROI,1);
% diffCoef = []; 
% dwellTotal = [];
% 
% for i = 1:CountROI
%     fileROI=sprintf(NameROI,i);
%     dataROI = csvread(fileROI,1,0);
% %     fileStat=sprintf(NameStat,i);
% %     statROI = csvread(fileStat,1,0);
%     % Filter trajectories by minimum number of steps and count length
%     dataROI(:,4) = dataROI(:,4) + 1; % Add one to count number of frames and not steps
%     dataROI=dataROI(any(dataROI(:,4)>=minStep,2),:);
%     index = dataROI(:,2);
%     dataFil=[];
%     for m = 1:size(index)
%         track = dataROI(any(dataROI(:,ci)==index(m),2),:);
%         dataFil = [dataFil ; track]; 
%     end
%     dwellTotal = [dwellTotal ; dataROI(:,4)]; % Merge all dwell times into a single file

%% Analysis density
   % density(i,1)=size(unique(dataFil(:,ci)),1)/area(i,1);
    
% %% Analysis of diffusion
%     dataFil(:,cx:cy) =  dataFil(:,cx:cy)/1000;
%     idx = unique(dataFil(:,ci));
%     diff = zeros(size(idx,1),1);
%     
%     for n = 1:size(idx,1)
%         p = find(dataFil(:,ci)==idx(n,1));
%         pi = p(1);
%         pf = p(end);
%         % Diffusion
%         msd = zeros(pf-pi,1);
%         for m = 1:pf-pi
%             dx = (dataFil(pi+m,cx) - dataFil(pi+m-1,cx))^2;
%             dy = (dataFil(pi+m,cy) - dataFil(pi+m-1,cy))^2;
%             msd(m,1) = dx + dy;
%         end
%         diff(n,1) = mean(msd)/(4*expT);
%     end
%     
%     % Merge results
%     diffCoef=[diffCoef ; diff];
%     
% 
% %% Analysis of PCA parameters
% % Parameter 1: Density
%      dataPCA(i,2) = density(i,1);
    
% % Parameter 2: Diffusion Coefficient
%      logDiff = log10(diff);
%      logDiff = logDiff(any(logDiff(:,1)~=-Inf,2),:); % Filterout struck particles with ratio of 0
%      pd = fitdist(logDiff,'Normal');
%      dataPCA(i,3) = pd.mu;
%      dataPCA(i,4) = pd.sigma;
    
% % Parameter 3: Total Distance
%     mDist = log10(statROI(:,3));
%     pd = fitdist(mDist,'Normal');
%     dataPCA(i,5) = pd.mu;
%     dataPCA(i,6) = pd.sigma;
% 
% % Parameter 4: Max distance
%     tDist = log10(statROI(:,4));
%     pd = fitdist(tDist,'Normal');
%     dataPCA(i,7) = pd.mu;
%     dataPCA(i,8) = pd.sigma;
%     
% % Parameter 5: Confinement Ratio
%     cRatio = log10(statROI(:,5));
%     cRatio = cRatio(any(cRatio(:,1)~=-Inf,2),:); % Filterout struck particles with ratio of 0
%     pd = fitdist(cRatio,'Normal');
%     dataPCA(i,9) = pd.mu;
%     dataPCA(i,10) = pd.sigma;
%     
% % Parameter 6: Mean Straigth Line Speed
%     mslSpeed = log10(statROI(:,6));
%     mslSpeed = mslSpeed(any(mslSpeed(:,1)~=-Inf,2),:);
%     pd = fitdist(mslSpeed,'Normal');
%     dataPCA(i,11) = pd.mu;
%     dataPCA(i,12) = pd.sigma;
% 
% % Parameter 7: Linearity of forward progression
%     lin = log10(statROI(:,7));
%     lin = lin(any(lin(:,1)~=-Inf,2),:);
%     pd = fitdist(lin,'Normal');
%     dataPCA(i,13) = pd.mu;
%     dataPCA(i,14) = pd.sigma;
% 
% % Parameter 8: Mean Directional Change Rate
%     mdcRate = statROI(:,8);
%     pd = fitdist(mdcRate,'Normal');
%     dataPCA(i,15) = pd.mu;
%     dataPCA(i,16) = pd.sigma;
%     
% % Parameter 9: Track duration/tau
%     h = transpose(histcounts(statROI(:,15), 'BinWidth', 1));
%     h(:,2) = minStep:size(h,1)+minStep-1;
%     h(:,3) = cumsum(h(:,1))/sum(h(:,1));
%     h=h(any(h(:,3)<=0.99,2),:);
%     
%     fitfun = fittype(@(a,t,x) a*exp(-x/t));
%     [fittedcurve, gof] = fit(h(:,2), h(:,1),fitfun, 'StartPoint', [h(1,2),h(1,1)]);
%     cv = coeffvalues(fittedcurve);
% 
%     dataPCA(i,17) = cv(1,2); % tau
% 
% % Parameter 10: Track displacement
%     tDisp = log10(statROI(:,18));
%     tDisp = tDisp(any(tDisp(:,1)~=-Inf,2),:);
%     pd = fitdist(tDisp,'Normal');
%     dataPCA(i,18) = pd.mu;
%     dataPCA(i,19) = pd.sigma;
%     
% % Parameter 11: Track mean speed
%     tmeanSpeed = log10(statROI(:,23));
%     tmeanSpeed = tmeanSpeed(any(tmeanSpeed(:,1)~=-Inf,2),:);
%     pd = fitdist(tmeanSpeed,'Normal');
%     dataPCA(i,20) = pd.mu;
%     dataPCA(i,21) = pd.sigma;
%     
% % Parameter 12: Track max speed
%     tmaxSpeed = log10(statROI(:,24));
%     tmaxSpeed = tmaxSpeed(any(tmaxSpeed(:,1)~=-Inf,2),:);
%     pd = fitdist(tmaxSpeed,'Normal');
%     dataPCA(i,22) = pd.mu;
%     dataPCA(i,23) = pd.sigma;
%     
% % Parameter 13: Track min speed 
%     tminSpeed = log10(statROI(:,25));
%     tminSpeed = tminSpeed(any(tminSpeed(:,1)~=-Inf,2),:);
%     pd = fitdist(tminSpeed,'Normal');
%     dataPCA(i,24) = pd.mu;
%     dataPCA(i,25) = pd.sigma;
%     
% % Parameter 14: Track median speed
%     tmedSpeed = log10(statROI(:,26));
%     tmedSpeed = tmedSpeed(any(tmedSpeed(:,1)~=-Inf,2),:);
%     pd = fitdist(tmedSpeed,'Normal');
%     dataPCA(i,26) = pd.mu;
%     dataPCA(i,27) = pd.sigma;
%     
% % Parameter 15: Track STD speed
%     tstdSpeed = log10(statROI(:,27));
%     tstdSpeed = tstdSpeed(any(tstdSpeed(:,1)~=-Inf,2),:);
%     pd = fitdist(tstdSpeed,'Normal');
%     dataPCA(i,28) = pd.mu;
%     dataPCA(i,29) = pd.sigma;
%  end
% 
% %% Save results PCA
%  csvwrite(namePCA,dataPCA);


%% Analysis of koff
%  h = transpose(histcounts(dwellTotal, 'BinWidth', 1));
%  h(:,2) = 1:size(h,1);
%  h = h(minStep:end,:);
%  h(:,3) = cumsum(h(:,1))/sum(h(:,1));
%  h=h(any(h(:,3)<=0.999,2),:);
%  
%  fitfun = fittype(@(a,t,x) a*exp(-x/t));
%  [fittedcurve, gof] = fit(h(:,2), h(:,1),fitfun, 'StartPoint', [h(1,2),h(1,1)]);
%  cv = coeffvalues(fittedcurve);
%  h(:,4) = cv(1,1)*exp(-h(:,2)/cv(1,2)) + 231.2*exp(-h(:,2)/15.2);
%  
% figure()
% bar(h(:,2),h(:,1)); hold on
% plot(h(:,2),h(:,4), 'LineWidth',2, 'Color', 'r');
% title('MDA231 B10')
% xlabel('Tracklength (frames)')
% ylabel('Frequency')
%  
% figure()
% scatter(h(:,2),log(h(:,1))); hold on
% plot(h(:,2),log(h(:,4)));
% axis([0 h(end,2) 0 (log(h(1,1))*1.2)])
% title('MDA231 B10')
% xlabel('Log Tracklength (frames)')
% ylabel('Frequency')
%  
%  koff = 1/(cv(1,2)*expT); %determining koff from the 1/tau times exposure time
%  koffError = koff*(gof.rmse/cv(1,1));
%  koffR = gof.adjrsquare;

%% Save results
 % csvwrite('Density_%d.csv',density);
%  csvwrite('DiffusionCoef.csv',diffCoef);
%  csvwrite('FitSecond.csv',h);