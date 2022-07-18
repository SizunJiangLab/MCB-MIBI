clear all
close all
tic;

%% INTERACTIVE: Set parameters
% Set 'paramsPath' to the location of the folder that contains the segmentationParams.mat file
paramsPath = '/Users/Desktop/folder/';

% Set 'clustersPath" to the location of the folder that contains the clusters.csv file
clustersPath = '/Users/Desktop/folder/';

% Set 'csv_file' to the name of the exported cluster file
csv_file = 'clusters.csv';

%% SCRIPT: cell coloring
%Load segmentationParams.mat and reads the segmentation mask information stored in newLmod) 
load([paramsPath,'segmentationParams.mat']);
labelNum = max(max(newLmod));
tempIm = zeros(size(newLmod,1),size(newLmod,2));

%Load cluster information
cluster_csv = readtable([clustersPath,csv_file]);
clusterID = cluster_csv.cluster;
clusterNum = max(max(clusterID));

%Make directory for saving cluster images
outputPath = [paramsPath, 'Cluster_output/'];
mkdir(outputPath);

%Modify newLmod such that each value is now the cluster value
%reshape(input, [(nrows * ncols)], 1]
reshaped_cellID = reshape(newLmod,[size(newLmod,1)*size(newLmod,2),1]); % reshapes newLmod matrix into a column vector

for k = 1:size(reshaped_cellID,1) %the (,1) informs the dimension of the reshape
    if reshaped_cellID(k) == 0
        clusterLmod_vector(k,1) = 0;
    else
        clusterLmod_vector(k,1) = clusterID(reshaped_cellID(k));
    end
end

clusterLmod = reshape(clusterLmod_vector,[size(newLmod,1),size(newLmod,2)]); % reshapes back into clusterLmod matrix
save([paramsPath,'/clusterParams.mat'],'clusterLmod')

%Make 1 color image per cluster
for k = 1:clusterNum
    for k1 = 1:size(clusterLmod_vector,1)
        if clusterLmod_vector(k1) == 0
            clusterLmod_channelvector{k}(k1,1) = 0;
        elseif clusterLmod_vector(k1) == k
            clusterLmod_channelvector{k}(k1,1) = 2;
        else
            clusterLmod_channelvector{k}(k1,1) = 1;
        end
    end
    clusterLmod_channel{k} = reshape(clusterLmod_channelvector{k},[size(newLmod,1),size(newLmod,2)]);
    clusterImage = imagesc(clusterLmod_channel{k});
    axis off;
    colormap bone;
    saveas(clusterImage,[outputPath,'colorMap_cluster',num2str(k),'.png']);
end

%Visualize all clusters
imagesc(clusterLmod)
colormap hot
toc;

%%%% The shortcut way
% for k = 1:clusterNum
%     a = clusterLmod;
%     a(a==0) = 0.5;
%     a(a>0.5 & a<k)=0;
%     a(a>0.5 & a>k)=0;
%     a(a==k)=2;
%     a(a==0)=1;
%     a(a==0.5)=0;
%     mapped_each{k} = a;
% end

%%%Make 1 color image per cluster
%%% This method does not work because I was looping over cellID, not cluster.
%%% Therefore, I needed another loop to loop over cell cluster.
%%% In the above code, k loops over cell clusters while k1 loops over cellID.
% for k = 1:size(reshaped_cellID,1)
%     if reshaped_cellID(k) == 0
%         clusterLmod_channelvector(k) = 0;
%     elseif reshaped_cellID(k) == k
%         clusterLmod_channelvector(k) = 2;
%     else
%         clusterLmod_channelvector(k) = 1;
%     end
%     clusterLmod_channel = reshape(clusterLmod_channelvector,[size(newLmod,1),size(newLmod,2)]);
%     clusterImage = imagesc(clusterLmod_channel);
%     imwrite(clusterImage,[outputPath,'colorMap_cluster',num2str(clusters),'.png']);
% end
