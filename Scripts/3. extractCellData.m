clear all
close all
tic;

%% INTERACTIVE: Set parameters
% Set 'mainPath' to the location of the folder that contains the raw data image files
mainPath = '/Users/Name/Folder/';

% Set 'mainFolder' to the name of the folder containing the raw data image files
mainFolder = 'FileName/';

% Indicate name of the DAPI channel to be used (without '.tif')
DAPI_name = 'Histone H3';

% If necessary, set 'fileSuffix' to the file format of the raw data image files (it is usually either .tif or .tiff).
fileSuffix = '.tiff';

% Set 'segMask' to the location of the segmentationMask.tif file
segMask = '/Users/Name/Folder/Cell_segmentation/';

%% SCRIPT: extract channel names
%Extract all the file names for all .tif files
Images = dir(strcat(mainPath,mainFolder,strcat('*',fileSuffix)));

for k = 1:size(Images,1)
    I_names{k,:} = Images(k).name;
end

%Obtain all the information from the file names
[pathstr,name,ext] = fileparts(I_names);

%Add header "Label" to the first row
header = {'Label'};
name = [header; name];

%Export channel names to CSV
writecell(name,[mainPath,'channelNames_analysis.csv']);

%% SCRIPT: extract single cell features and save to .FCS and .csv
%Read channelNames_analysis.csv file and obtain the channel names
channelNames = dataset('File',[mainPath,'channelNames_analysis.csv',],'Delimiter',',');
clusterChannels = channelNames.Label;
[~, clusterChannelsInds] = ismember(clusterChannels,channelNames.Label);

for p=1
    disp(['p',num2str(p)]);
    pointNumber = p;
    % load tiffs to recreate countsNoNoise
    for i = 1:length(channelNames.Label)
        t = imread([mainPath,mainFolder, channelNames.Label{i}, fileSuffix]); 
        d = double(t);
        countsNoNoise(:,:,i) = d;
    end

    %Load segmentation results from segmentationMask.tif
    newLmod = imread([segMask,'segmentationMask.tif']);
    labelNum = max(max(newLmod));
    labelIdentityNew = ones(1,labelNum);
    channelNum = length(channelNames);

    %Reshape countsNoNoise into #pixels x #images matrix
    countsReshape= reshape(countsNoNoise,size(countsNoNoise,1)*size(countsNoNoise,2),channelNum);
    
    %Segment cells by regionprops and reshapes 
    stats = regionprops(newLmod,'Area','PixelIdxList');

    %Generate (x,y) coordinates of each cell
    coordinates = regionprops(newLmod,'Centroid')

    %Make a data matrix comprising #labels x #markers (including one more marker for cell size)
    data = zeros(labelNum,channelNum);
    dataScaleSize = zeros(labelNum,channelNum);
    cellSizes = zeros(labelNum,1);
    
    %For each label, extract the information from the countsReshape matrix
    for i=1:labelNum
        currData = countsReshape(stats(i).PixelIdxList,:);
        data(i,1:channelNum) = sum(currData,1);
        dataScaleSize(i,1:channelNum) = sum(currData,1) / stats(i).Area;
        cellSizes(i) = stats(i).Area;
    end
       
    %Get the final information only for the labels with:
    % 1. Positive nuclear identity (indicates cells)
    % 2. Enough information in the clustering channels to be clustered
    sumDataScaleSizeInClusterChannels = sum(dataScaleSize(:,clusterChannelsInds),2);
    labelIdentityNew(sumDataScaleSizeInClusterChannels<0.1) = 2; % this line may be troublesome as the 0.1 value is for MIBI    
    labelVec=find(labelIdentityNew==1)';

    %Get the cell sizes from the selected cells
    cellSizesVec = cellSizes(labelIdentityNew==1);

    %Store the features from all cells into the dataL vector
    dataCells = data(labelIdentityNew==1,:);
    dataL = [labelVec,cellSizesVec,dataCells];

    %Store the features from all cells (but scaled to cell size) into the dataScaleSizeL vector
    dataScaleSizeCells = dataScaleSize(labelIdentityNew==1,:);    
    dataScaleSizeL = [labelVec,cellSizesVec,dataScaleSizeCells];
    
    %% Save single cell features
    % Feel free to change the 'pathResuls' file name, where the all the extracted cell data will be saved
    pathResults = ['FCS_output'];    
    outputPath = [mainPath, '/', pathResults];
    mkdir(outputPath);

    %Channel labels for the FCS output
    channelLabelsForFCS = ['cellLabelInImage';'cellSize';channelNames.Label;'PointNum'];

    %Save fcs
    TEXT.PnS = channelLabelsForFCS;
    TEXT.PnN = channelLabelsForFCS;
    save([mainPath,'/segmentationParams.mat'],'newLmod','labelIdentityNew');
    save([outputPath,'/cellData.mat'],'labelIdentityNew','labelVec','cellSizesVec','dataCells','dataScaleSizeCells','channelLabelsForFCS'); 
    
    %Feel free to change the file names between the quotation marks
    writeFCS([outputPath,'/data.fcs'],dataL,TEXT);
    writeFCS([outputPath,'/dataScaleSize.fcs'],dataScaleSizeL,TEXT);
    csvwrite([outputPath,'/data.csv'],dataL);
    csvwrite([outputPath,'/dataScaleSize.csv'],dataScaleSizeL); 
end
toc;
