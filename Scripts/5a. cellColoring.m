clear all
close all
tic;

%% INTERACTIVE: Set parameters
% Set 'mainPath' to the location of the folder that contains the segmentationParams.mat file
paramsPath = '/Users/Desktop/folder/';

%% Define colors
% This is the 20 color codes the script will be using for cell coloring. Feel free to change the colors if desired.
C ={'177FFA','55B73C','FFB921','F2F200','97AA37','2CF40E','00ECEC', 'E20550','7300CC', 'AE2CCA', ...
    '81FFE0','D7FF9F','F9FAFC','C7A1FF', 'CFCAAF','ABBBD3','FFBBFF','2400FC','C9E8DE','B1C0FF'};

%% SCRIPT: cell coloring
%Load segmentationParams.mat and reads the segmentation mask information stored in newLmod) 
load([paramsPath,'segmentationParams.mat']);
labelNum = max(max(newLmod));
tempIm = zeros(size(newLmod,1),size(newLmod,2));

%Obtains the colors to be used
colorHex = string(C);
color = zeros(20,3);
for i = 1:20
    color(i,:) = hex2rgb(colorHex{i});
end

%Colors all the cells. This takes a while to run.
for cells = 1:labelNum
    if labelIdentityNew(cells) == 1 %Stains cells with color
        tempIm = colorTheMap(cells,tempIm,newLmod,color(randi([1 20]),:));   
    elseif labelIdentityNew(cells) ~= 1 %Stains non-cells as black
        tempIm = colorTheMap(cells,tempIm,newLmod,[0 0 0]);
    end
end

%Display coloring map and save coloring map
figure;
imshow(tempIm);
imwrite(tempIm,[paramsPath,'colorMap_randomColor.png']);
toc;

%% Define functions for this script
%Assist function for map coloring
function tempIm = colorTheMap(cell,tempIm,newLmod,colorSet)
    [tempRow,tempCol] = find(newLmod==cell);
    for j = 1:length(tempRow)
        tempIm(tempRow(j),tempCol(j),1) = colorSet(1);
        tempIm(tempRow(j),tempCol(j),2) = colorSet(2);
        tempIm(tempRow(j),tempCol(j),3) = colorSet(3);
    end
end

%Assist function for color vector
function rgb = hex2rgb(hexString)
    r = double(hex2dec(hexString(1:2)))/255;
    g = double(hex2dec(hexString(3:4)))/255;
    b = double(hex2dec(hexString(5:6)))/255;
	rgb = [r, g, b];
end
