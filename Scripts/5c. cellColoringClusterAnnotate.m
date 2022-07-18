clear all
close all
tic;

%% INTERACTIVE: Set parameters
% Set 'paramsPath' to the location of the folder that contains the clusterParams.mat file
paramsPath = '/Users/Desktop/folder/';

% Set 'annotatePath" to the location of the folder that contains the annotated.csv file
annotatePath = '/Users/Desktop/folder/';

% Set 'csv_file' to the name of the annotated file
csv_file = 'annotated.csv';

% Scroll to the bottom of the script to modify image parameters

%% SCRIPT: cell coloring
%Load clusterParams.mat and reads the segmentation mask information stored in newLmod) 
load([paramsPath,'clusterParams.mat']);
labelNum = max(max(clusterLmod));
tempIm = zeros(size(clusterLmod,1),size(clusterLmod,2));

%Load annotation information
annotate_csv = readtable([annotatePath,csv_file]);
annotateID = annotate_csv.celltype;
annotateNum = max(max(annotateID));

%Make directory for saving annotated cluster image
outputPath = [paramsPath, 'Annotate_output/'];
mkdir(outputPath);

%Modify newLmod such that each value is now the cluster value
%reshape(input, [(nrows * ncols)], 1]
reshaped_annotateID = reshape(clusterLmod,[size(clusterLmod,1)*size(clusterLmod,2),1]); % reshapes clusterLmod matrix into a column vector

for k = 1:size(reshaped_annotateID,1) %the (,1) informs the dimension of the reshape
    if reshaped_annotateID(k) == 0
        annotateLmod_vector(k,1) = 0;
    else
        annotateLmod_vector(k,1) = annotateID(reshaped_annotateID(k));
    end
end

annotateLmod = reshape(annotateLmod_vector,[size(clusterLmod,1),size(clusterLmod,2)]); % reshapes back into annotateLmod matrix
save([paramsPath,'/annotateParams.mat'],'annotateLmod')

%Make 1 color image per cluster
for k = 1:annotateNum
    for k1 = 1:size(annotateLmod_vector,1)
        if annotateLmod_vector(k1) == 0
            annotateLmod_channelvector{k}(k1,1) = 0;
        elseif annotateLmod_vector(k1) == k
            annotateLmod_channelvector{k}(k1,1) = 2;
        else
            annotateLmod_channelvector{k}(k1,1) = 1;
        end
    end
    annotateLmod_channel{k} = reshape(annotateLmod_channelvector{k},[size(clusterLmod,1),size(clusterLmod,2)]);
    annotateImage = imagesc(annotateLmod_channel{k});
    axis off;
    colormap bone;
    saveas(annotateImage,[outputPath,'colorMap_annotate',num2str(k),'.png']);
end

%% INTERACTIVE: Export image
% These are the color codes the script will be using for cell coloring. Feel free to change the colors if desired.
% Tonsil: 8 cell types = 8 colors
Colors = {'4DAF4A','E41A1C','377EB8','984EA3','FF7F00','FFFF33','A65628','F781BF'};
%%% Additional random colors: {'0B6DE0','8E7CC3','FFB921','FFFF0B','97AA37','2CF40E','00ECEC', 'E20550'}

% Indicate the name of all cell types
CellName = ["B cell", "CD8 T cell", "CD4 T cell", "Treg", "NK cell", "M1 MΦ", "M2 MΦ", "Dendritic cell"];

%Obtains the colors to be used
C = ['000000', Colors]; % adding black (#000000) to color background
colorHex = string(C);
color = zeros(length(C),3);
for i = 1:length(C)
    color(i,:) = hex2rgb(colorHex{i});
end

%Visualize all clusters, and save image to outputPath
colorImage = imagesc(annotateLmod);
colormap(color);
axis off;
saveas(colorImage,[outputPath,'colorImage.png'])

%Visualize all clusters and save image with labels to outputPath
CellName = ["", CellName]; % adding empty ("") to label background
colorImageLabel = imagesc(annotateLmod);
colormap(color);
colorbar( ...
    'Ticks',0:length(C), ...
    'TickLabels',CellName);
saveas(colorImageLabel,[outputPath,'colorImageLabel.png'])
toc;

%% Define functions for this script
%Assist function for color vector
function rgb = hex2rgb(hexString)
    r = double(hex2dec(hexString(1:2)))/255;
    g = double(hex2dec(hexString(3:4)))/255;
    b = double(hex2dec(hexString(5:6)))/255;
	rgb = [r, g, b];
end
