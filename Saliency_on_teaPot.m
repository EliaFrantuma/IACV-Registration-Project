%% Start - Up
clear;
close all;
%% Load + normals

layers = 3;



for i = 1:12
    readptCloud(i) = pcdenoise( pcread(strcat(strcat("Armadillo_scans\ArmadilloStand_",string((i-1)*30)),".ply")) );
    readptCloud(i) = pointCloud(unique( readptCloud(i).Location, 'rows'));
    readptCloud(i) = pcdownsample(readptCloud(i),'gridAverage',0.001);

end

ptCloud(1) = readptCloud(4);
ptCloud(2) = readptCloud(5);


maxDistance = norm([ptCloud(1).XLimits(1) ptCloud(1).YLimits(1) ptCloud(1).ZLimits(1)] - [ptCloud(1).XLimits(2) ptCloud(1).YLimits(2) ptCloud(1).ZLimits(2)]);

baseSigma = maxDistance/30;


normals = pcnormals(ptCloud(1),100);
normals2 = pcnormals(ptCloud(2),100);


wrongIndices1 = cleanCloud(ptCloud(1),baseSigma);
wrongIndices2 = cleanCloud(ptCloud(2),baseSigma);

%      figure
%      pcshow(ptCloudNoBorder(1))

%% 



% for c = 1:layers
%     %smoothedClouds(c) = pointCloud(imgaussfilt3(ptCloud(1).Location,(2^c)/denom));
%     for point = 1:ptCloud(1).Count
%         assert(size(findNeighborsInRadius(ptCloud(1),ptCloud(1).Location(point,:),(2^(c-1))*baseSigma),1) == size(findNeighborsInRadius(ptCloud(2),ptCloud(2).Location(point,:),(2^(c-1))*baseSigma),1));
%     end
% end

%% Gaussian 3D filtering

smoothedClouds = pointCloud.empty(0,layers);
smoothedClouds2 = pointCloud.empty(0,layers);

smoothedLocation = zeros(ptCloud(1).Count,3);
smoothedLocation2 = zeros(ptCloud(2).Count,3);

smoothedClouds(1) = ptCloud(1);
smoothedClouds2(1) = ptCloud(2);


for c = 2:layers
    %smoothedClouds(c) = pointCloud(imgaussfilt3(ptCloud(1).Location,(2^c)/denom));
    for point = 1:ptCloud(1).Count
        smoothedLocation(point,:) = smoothPointInCloud(ptCloud(1),ptCloud(1).Location(point,:),((c-1))*baseSigma);
    end
    smoothedClouds(c) = pointCloud(smoothedLocation);
    for point = 1:ptCloud(2).Count
        smoothedLocation2(point,:) = smoothPointInCloud(ptCloud(2),ptCloud(2).Location(point,:),((c-1))*baseSigma);
    end
    smoothedClouds2(c) = pointCloud(smoothedLocation2);
end


%% Gaussian Projection

gaussProj = pointCloud.empty(0,layers);
gaussProj2 = pointCloud.empty(0,layers);
gaussNormals = zeros(layers, ptCloud(1).Count, 3);
gaussNormals2 = zeros(layers, ptCloud(2).Count, 3);

for c = 1:layers
    temploc = zeros(ptCloud(1).Count,3);
    temploc2 = zeros(ptCloud(2).Count,3);

    for i = 1 : size(ptCloud(1).Location,1)
        temploc(i,:) = ptCloud(1).Location(i,:) + dot((smoothedClouds(c).Location(i,:)-ptCloud(1).Location(i,:)),normals(i,:))*normals(i,:);
    end
    for i = 1 : size(ptCloud(2).Location,1)
        temploc2(i,:) = ptCloud(2).Location(i,:) + dot((smoothedClouds2(c).Location(i,:)-ptCloud(2).Location(i,:)),normals2(i,:))*normals2(i,:);
    end
    gaussProj(c) = pointCloud(temploc);
    gaussProj2(c) = pointCloud(temploc2);
    
    gaussNormals(c,:,:) = pcnormals(gaussProj(c),12);
    gaussNormals2(c,:,:) = pcnormals(gaussProj2(c),12);
end

%% Saliency maps

for c = 1:layers-1
    saliencyMaps(c).mat = gaussProj(c).Location - gaussProj(c+1).Location;
    normMat(c).mat = vecnorm((saliencyMaps(c).mat)')';
%     normalsMat = reshape(gaussNormals(c,:,:),size(gaussNormals(c,:,:),2),size(gaussNormals(c,:,:),3)) * reshape(gaussNormals(c+1,:,:),size(gaussNormals(c+1,:,:),2),size(gaussNormals(c+1,:,:),3))';
%     normMat(c).mat = normMat(c).mat .* normalsMat(); 
end

for c = 1:layers-1
    saliencyMaps2(c).mat = gaussProj2(c).Location - gaussProj2(c+1).Location;
    normMat2(c).mat = vecnorm((saliencyMaps2(c).mat)')';
end

%% Maximum Extraction 1

count=1;

figure
pcshow(ptCloud(1).Location,normMat(1).mat)
hold on

for c = 1:layers-1
for point = 1:size(ptCloud(1).Location,1)
    if ~ismember(point,wrongIndices1)
    [indices,dists] = findNeighborsInRadius(ptCloud(1),ptCloud(1).Location(point,:),((c))*baseSigma,'Sort',true);
    indices = setdiff(indices,wrongIndices1,'stable');
    flag = true;
    for i = 2:size(indices)
        if normMat(c).mat(indices(i)) >= normMat(c).mat(point)
            flag = false;
        end
    end
    if flag %&& size(indices,1) > 60*((c))
        maxIndices(count,c) = point;
        plot3(ptCloud(1).Location(point,1),ptCloud(1).Location(point,2),ptCloud(1).Location(point,3),'*r')
        count=count+1;
    end
    end
end
count=1;
end

plot3(ptCloud(1).Location(1,1),ptCloud(1).Location(1,2),ptCloud(1).Location(1,3),'*g')
hold off

%% Maximum extraction 2


count=1;

figure
pcshow(ptCloud(2).Location,normMat2(1).mat)
hold on

for c = 1:layers-1
for point = 1:size(ptCloud(2).Location,1)
    if ~ismember(point,wrongIndices2)
    [indices,dists] = findNeighborsInRadius(ptCloud(2),ptCloud(2).Location(point,:),((c))*baseSigma,'Sort',true);
    indices = setdiff(indices,wrongIndices2,'stable');
    flag = true;
    for i = 2:size(indices)
        if normMat2(c).mat(indices(i)) >= normMat2(c).mat(point)
           flag = false;
        end
    end
    if flag% && size(indices,1) > 60*(c)
        maxIndices2(count,c) = point;
        plot3(ptCloud(2).Location(point,1),ptCloud(2).Location(point,2),ptCloud(2).Location(point,3),'*r')
        count=count+1;
    end
    end
end
count=1;
end

plot3(ptCloud(2).Location(1,1),ptCloud(2).Location(1,2),ptCloud(2).Location(1,3),'*g')

hold off


%% Descriptor generation

M=4;
L=8;

saliencyNorms = zeros(size(normMat(1).mat,1),layers-1);

for c =1:layers-1
    for i = 1:size(normMat(1).mat,1)
        saliencyNorms(i,c)=normMat(c).mat(i);
    end
end

maxIndicesVector = reshape(maxIndices,[],1);
descriptorsVector=[]; 

missing = 0;

for c = 1:layers-1

for index=1:size(maxIndices,1)
    if maxIndices(index,c)>0
        descriptorsVector(index,c).zDir = normals(maxIndices(index,c),:);
        xyMat = null(descriptorsVector(index,c).zDir(:).');
        if abs(cross(xyMat(:,1),xyMat(:,2)) - descriptorsVector(index,c).zDir(:) ) < 0.0001
            descriptorsVector(index,c).xDir = xyMat(:,1);
            descriptorsVector(index,c).yDir = xyMat(:,2);
        else
            descriptorsVector(index,c).xDir = xyMat(:,2);
            descriptorsVector(index,c).yDir = xyMat(:,1);
        end

        [indicesOfNeigh,distances] = findNeighborsInRadius(ptCloud(1),ptCloud(1).Location(maxIndices(index,c),:),((c+1))*baseSigma,'Sort',true);
        neighborsGrid = zeros(size(indicesOfNeigh,1),3); %indices legend : M=1, L=2, CloudIndex=3
        for neighIndex = 1:size(indicesOfNeigh)
            diffVector = ptCloud(1).Location(maxIndices(index,c),:) - ptCloud(1).Location(indicesOfNeigh(neighIndex),:);      
            sigma = ((c+1))*baseSigma;
            neighborsGrid(neighIndex,1) = ceil((vecnorm(diffVector)*M)/sigma);
            phi = acos(dot(diffVector,descriptorsVector(index,c).zDir));
            normalOnPlane = (diffVector-(vecnorm(diffVector)*cos(phi) * descriptorsVector(index,c).zDir))/vecnorm(diffVector-(vecnorm(diffVector)*cos(phi) * descriptorsVector(index,c).zDir));
            normalsAndY=dot(normalOnPlane,descriptorsVector(index,c).yDir);
            if normalsAndY >= 0
                theta = acos(dot(normalOnPlane,descriptorsVector(index,c).xDir));
            else
                theta = 2*pi - acos(dot(normalOnPlane,descriptorsVector(index,c).xDir));
            end
            
            
            neighborsGrid(neighIndex,2) = ceil((theta*L)/(2*pi));
            neighborsGrid(neighIndex,3) = indicesOfNeigh(neighIndex);
        end
        neighborsGrid(neighborsGrid(:,1) > M,1) = M;
        neighborsGrid(neighborsGrid(:,1) < 1,1) = 1;
        neighborsGrid(isnan(neighborsGrid(:,2)),2) = 1;
        neighborsGrid = real(neighborsGrid);

        singleDescriptorGrid = zeros(M,L,2);
        singleDescriptorGridMask = ones(M,L,2);
        
        for range = 1:M
            for angle = 1:L
                %normalMean = mean(normals(M(:,1) == someValue, 3)); 
                  allM = neighborsGrid(:,1); % comma separated list expansion 
                  tfM = allM == range;
                  indexM = find(tfM) ;
                  allL = neighborsGrid(:,2); % comma separated list expansion 
                  tfL = allL == angle;
                  indexL = find(tfL) ;
                  
                  pointsInCell = intersect(indexM,indexL);
                  
                  averageNorm = [0,0,0];
                  averageSaliency = 0;
                  
                  if size(pointsInCell)>0
                    for i = 1:size(pointsInCell)
                        averageNorm = averageNorm + normals(neighborsGrid(pointsInCell(i),3),:);
                        averageSaliency = averageSaliency + saliencyNorms(neighborsGrid(pointsInCell(i),3),c);
                    end
                    averageNorm = averageNorm / norm(averageNorm);
                    averageSaliency = averageSaliency / size(pointsInCell,1);
                  else
                      singleDescriptorGridMask(range,angle,:) = [0 0];
                    
                  end
                  
                  
                  deltaNorm = 1.0 - abs(dot(averageNorm,normals(maxIndices(index,c),:)));
                  deltaSaliency = 1.0 - averageSaliency/saliencyNorms(maxIndices(index,c),c);
                  singleDescriptorGrid(range,angle,:) = [deltaNorm, deltaSaliency];
            end
        end
        
       % if maxIndicesVector(index)==23748 || maxIndicesVector(index)==25067 || maxIndicesVector(index)==5506
        maxIndicesVectorDescriptor(index-missing,c).descriptor = singleDescriptorGrid;
        maxIndicesVectorDescriptor(index-missing,c).mask = singleDescriptorGridMask;
        maxIndicesVectorDescriptor(index-missing,c).index = maxIndices(index,c);
       % else
      %      missing = missing+1;
      %  end
        
    end
end
end

%% Descriptor generation 2

saliencyNorms2 = zeros(size(saliencyMaps2(1).mat,1),layers-1);

for c =1:layers-1
    for i = 1:size(saliencyMaps2(1).mat,1)
        saliencyNorms2(i,c)=normMat2(c).mat(i);
    end
end

maxIndicesVector2 = reshape(maxIndices2,[],1);
descriptorsVector2=[]; 

missing = 0;
for c = 1:layers-1

for index=1:size(maxIndices2)
    if maxIndices2(index,c)>0
        descriptorsVector2(index,c).zDir = normals2(maxIndices2(index,c),:);
        xyMat = null(descriptorsVector2(index,c).zDir(:).');
        if abs(cross(xyMat(:,1),xyMat(:,2)) - descriptorsVector2(index,c).zDir(:) ) < 0.0001
            descriptorsVector2(index,c).xDir = xyMat(:,1);
            descriptorsVector2(index,c).yDir = xyMat(:,2);
        else
            descriptorsVector2(index,c).xDir = xyMat(:,2);
            descriptorsVector2(index,c).yDir = xyMat(:,1);
        end
        [indicesOfNeigh,distances] = findNeighborsInRadius(ptCloud(2),ptCloud(2).Location(maxIndices2(index,c),:),((c+1))*baseSigma,'Sort',true);
        neighborsGrid = zeros(size(indicesOfNeigh,1),3); %indices legend : M=1, L=2, CloudIndex=3
        for neighIndex = 1:size(indicesOfNeigh)
            diffVector = ptCloud(2).Location(maxIndices2(index,c),:) - ptCloud(2).Location(indicesOfNeigh(neighIndex),:);      
            sigma = ((c+1))*baseSigma;
            neighborsGrid(neighIndex,1) = ceil((vecnorm(diffVector)*M)/sigma);

            phi = acos(dot(diffVector,descriptorsVector2(index,c).zDir));
            normalOnPlane = (diffVector-(vecnorm(diffVector)*cos(phi) * descriptorsVector2(index,c).zDir))/vecnorm(diffVector-(vecnorm(diffVector)*cos(phi) * descriptorsVector2(index,c).zDir));
            normalsAndY=dot(normalOnPlane,descriptorsVector2(index,c).yDir);
            if normalsAndY >= 0
                theta = acos(dot(normalOnPlane,descriptorsVector2(index,c).xDir));
            else
                theta = 2*pi - acos(dot(normalOnPlane,descriptorsVector2(index,c).xDir));
            end
            
            neighborsGrid(neighIndex,2) = ceil((theta*L)/(2*pi));
            neighborsGrid(neighIndex,3) = indicesOfNeigh(neighIndex);
        end
        
        neighborsGrid(neighborsGrid(:,1) > M,1) = M;
        neighborsGrid(neighborsGrid(:,1) < 1,1) = 1;
        neighborsGrid(isnan(neighborsGrid(:,2)),2) = 1;
        neighborsGrid = real(neighborsGrid);
        
        singleDescriptorGrid2 = zeros(M,L,2);
        singleDescriptorGridMask2 = ones(M,L,2);

        
        for range = 1:M
            for angle = 1:L
                %normalMean = mean(normals(M(:,1) == someValue, 3)); 
                  allM = neighborsGrid(:,1); % comma separated list expansion 
                  tfM = allM == range;
                  indexM = find(tfM) ;
                  allL = neighborsGrid(:,2); % comma separated list expansion 
                  tfL = allL == angle;
                  indexL = find(tfL) ;
                  
                  pointsInCell = intersect(indexM,indexL);
                  
                  averageNorm = [0,0,0];
                  averageSaliency = 0;
                  
                  if size(pointsInCell)>0
                    for i = 1:size(pointsInCell)
                        averageNorm = averageNorm + normals2(neighborsGrid(pointsInCell(i),3),:);
                        averageSaliency = averageSaliency + saliencyNorms2(neighborsGrid(pointsInCell(i),3),c);
                    end
                    averageNorm = averageNorm / norm(averageNorm);
                    averageSaliency = averageSaliency / size(pointsInCell,1);
                  else
                      singleDescriptorGridMask2(range,angle,:) = [0 0];
                    
                  end
                  
                  deltaNorm = 1.0 - abs(dot(averageNorm,normals2(maxIndices2(index,c),:)));
                  deltaSaliency = 1.0 - averageSaliency/ saliencyNorms2(maxIndices2(index,c),c);
                  singleDescriptorGrid2(range,angle,:) = [deltaNorm, deltaSaliency];
            end
        end
        
       % if maxIndicesVector2(index)==10572 || maxIndicesVector2(index)==11891 || maxIndicesVector2(index)==2806
        maxIndicesVectorDescriptor2(index-missing,c).descriptor = singleDescriptorGrid2;
        maxIndicesVectorDescriptor2(index-missing,c).mask = singleDescriptorGridMask2;

        maxIndicesVectorDescriptor2(index-missing,c).index = maxIndices2(index,c);
       % else
       %     missing = missing +1;
       % end
        
    end
end
end


%% Feature Matching

length1=size(maxIndices,1);
length2=size(maxIndices2,1);
maxCompleteScore = zeros(size(maxIndicesVectorDescriptor,2),size(maxIndicesVectorDescriptor2,2));
maxCompleteScore = zeros(length1,length2,layers-1);

for c = 1:layers-1
for s = 1:length1
    for d = 1:length2
        maxCompleteScore(s,d,c) = 0;            
        if ~isempty(maxIndicesVectorDescriptor2(d,c).mask(:,:,1)) && ~isempty(maxIndicesVectorDescriptor(s,c).mask(:,:,1))
        for l_bar = 1:L
            normalScore = maxIndicesVectorDescriptor(s,c).mask(:,:,1).*circshift(maxIndicesVectorDescriptor2(d,c).mask(:,:,1),l_bar,2).*(1 - abs(maxIndicesVectorDescriptor(s,c).descriptor(:,:,1) - circshift(maxIndicesVectorDescriptor2(d,c).descriptor(:,:,1),l_bar,2)));
            saliencyScore = maxIndicesVectorDescriptor(s,c).mask(:,:,2).*circshift(maxIndicesVectorDescriptor2(d,c).mask(:,:,2),l_bar,2).*(1 - abs(maxIndicesVectorDescriptor(s,c).descriptor(:,:,2) - circshift(maxIndicesVectorDescriptor2(d,c).descriptor(:,:,2),l_bar,2)));
            
            normalScoreArray = normalScore';
            normalScoreArray = normalScoreArray(:)';
            
            saliencyScoreArray = saliencyScore';
            saliencyScoreArray = saliencyScoreArray(:)';
            
            tempCompleteSD = dot(normalScoreArray,saliencyScoreArray);
            if tempCompleteSD > maxCompleteScore(s,d,c)
                maxCompleteScore(s,d,c) = tempCompleteSD;
            end
        end
        end
        
    end
end
end

%% Selection of top 5 features

%matchedFeatures = zeros(100,2);

for index=1:5
    [mxv,idx] = max(maxCompleteScore(:));
    [row,col,layer] = ind2sub(size(maxCompleteScore),idx);
    
    row = row(1);
    col = col(1);
    layer = layer(1);
    matchedFeatures(index,:) = [row col layer];
    maxCompleteScore(row,:,layer) = 0;
    maxCompleteScore(:,col,layer) = 0;
    if all(all(ismember(maxCompleteScore, max(maxCompleteScore(:))))) 
        break
    end

end

%% Get Point clouds from matched points

firstPoints = matchedFeatures(:,1) + (matchedFeatures(:,3)-1)*length1;
secondPoints = matchedFeatures(:,2) + (matchedFeatures(:,3)-1)*length2;


fixedCloud = ptCloud(1).Location(maxIndicesVector(firstPoints),:);
changingCloud = ptCloud(2).Location(maxIndicesVector2(secondPoints),:);



%% RANSAC

[tformEst,inlierIndex] = estimateGeometricTransform3D(changingCloud,fixedCloud,'rigid','maxNumTrials',100000,'Confidence',99); 

inliersPtCloud2 = pctransform(ptCloud(2),tformEst).Location;
inliersPtCloud1 = ptCloud(1).Location;

% inliersPtCloud2 = pctransform(pointCloud(changingCloud),tformEst).Location;
% inliersPtCloud1 = fixedCloud;


figure
firstPtCloud = pointCloud(inliersPtCloud1);
secondPtCloud = pointCloud(inliersPtCloud2);
pcshowpair(firstPtCloud,secondPtCloud)
title('Aligned point clouds')

% figure
% hold on
% firstPtCloud1 = ptCloud(1);
% pcshow(firstPtCloud1)
% plot3(fixedCloud(1,1),fixedCloud(1,2),fixedCloud(1,3),'*b')
% title('Aligned point clouds')
% hold off
% 
% figure
% hold on
% firstPtCloud1 = ptCloud(2);
% pcshow(firstPtCloud1)
% plot3(changingCloud(1,1),changingCloud(1,2),changingCloud(1,3),'*b')
% title('Aligned point clouds')
% hold off



%% ICP

ptCloudRef = firstPtCloud;
ptCloudCurrent = secondPtCloud;

gridSize = 0.0001;
fixed = pcdownsample(ptCloudRef, 'gridAverage', gridSize);
moving = pcdownsample(ptCloudCurrent, 'gridAverage', gridSize);

tform = pcregistericp(moving, fixed, 'Metric','pointToPlane','Extrapolate', true);
ptCloudAligned = pctransform(ptCloudCurrent,tform);

mergeSize = 0.0015;
ptCloudScene = pcmerge(ptCloudRef, ptCloudAligned, mergeSize);

figure
hAxes = pcshow(ptCloudScene, 'VerticalAxis','Y', 'VerticalAxisDir', 'Down');
title('Updated world scene')
% Set the axes property for faster rendering
hAxes.CameraViewAngleMode = 'auto';
hScatter = hAxes.Children;

t12 = tformEst.T * tform.T;
% Store the transformation object that accumulates the transformation.
accumTform = rigid3d(t12);

for iteratore = 1:10

iterator = mod(iteratore+5,13);
if iteratore >= 8
    iterator = iterator+1;
end
if iterator == 1
    ptCloud(1) = readptCloud(12);
    ptCloud(2) = readptCloud(1);
else
    ptCloud(1) = readptCloud(iterator-1);
    ptCloud(2) = readptCloud(iterator);
end

maxIndices = [];
maxIndices2 = [];
maxIndicesVectorDescriptor = [];
maxIndicesVectorDescriptor2 = [];

if iterator == 7
    denom = 25;
elseif iterator == 9
    denom = 20;
else
    denom = 30;
end






maxDistance = norm([ptCloud(1).XLimits(1) ptCloud(1).YLimits(1) ptCloud(1).ZLimits(1)] - [ptCloud(1).XLimits(2) ptCloud(1).YLimits(2) ptCloud(1).ZLimits(2)]);

baseSigma = maxDistance/denom;


normals = pcnormals(ptCloud(1),100);
normals2 = pcnormals(ptCloud(2),100);


wrongIndices1 = cleanCloud(ptCloud(1),baseSigma);
wrongIndices2 = cleanCloud(ptCloud(2),baseSigma);
%% Gaussian 3D filtering

smoothedClouds = pointCloud.empty(0,layers);
smoothedClouds2 = pointCloud.empty(0,layers);

smoothedLocation = zeros(ptCloud(1).Count,3);
smoothedLocation2 = zeros(ptCloud(2).Count,3);

smoothedClouds(1) = ptCloud(1);
smoothedClouds2(1) = ptCloud(2);


for c = 2:layers
    %smoothedClouds(c) = pointCloud(imgaussfilt3(ptCloud(1).Location,(2^c)/denom));
    for point = 1:ptCloud(1).Count
        smoothedLocation(point,:) = smoothPointInCloud(ptCloud(1),ptCloud(1).Location(point,:),((c-1))*baseSigma);
    end
    smoothedClouds(c) = pointCloud(smoothedLocation);
    for point = 1:ptCloud(2).Count
        smoothedLocation2(point,:) = smoothPointInCloud(ptCloud(2),ptCloud(2).Location(point,:),((c-1))*baseSigma);
    end
    smoothedClouds2(c) = pointCloud(smoothedLocation2);
end


%% Gaussian Projection

gaussProj = pointCloud.empty(0,layers);
gaussProj2 = pointCloud.empty(0,layers);
gaussNormals = zeros(layers, ptCloud(1).Count, 3);
gaussNormals2 = zeros(layers, ptCloud(2).Count, 3);

for c = 1:layers
    temploc = zeros(ptCloud(1).Count,3);
    temploc2 = zeros(ptCloud(2).Count,3);

    for i = 1 : size(ptCloud(1).Location,1)
        temploc(i,:) = ptCloud(1).Location(i,:) + dot((smoothedClouds(c).Location(i,:)-ptCloud(1).Location(i,:)),normals(i,:))*normals(i,:);
    end
    for i = 1 : size(ptCloud(2).Location,1)
        temploc2(i,:) = ptCloud(2).Location(i,:) + dot((smoothedClouds2(c).Location(i,:)-ptCloud(2).Location(i,:)),normals2(i,:))*normals2(i,:);
    end
    gaussProj(c) = pointCloud(temploc);
    gaussProj2(c) = pointCloud(temploc2);
    
    gaussNormals(c,:,:) = pcnormals(gaussProj(c),12);
    gaussNormals2(c,:,:) = pcnormals(gaussProj2(c),12);
end

%% Saliency maps

for c = 1:layers-1
    saliencyMaps(c).mat = gaussProj(c).Location - gaussProj(c+1).Location;
    normMat(c).mat = vecnorm((saliencyMaps(c).mat)')';
%     normalsMat = reshape(gaussNormals(c,:,:),size(gaussNormals(c,:,:),2),size(gaussNormals(c,:,:),3)) * reshape(gaussNormals(c+1,:,:),size(gaussNormals(c+1,:,:),2),size(gaussNormals(c+1,:,:),3))';
%     normMat(c).mat = normMat(c).mat .* normalsMat(); 
end

for c = 1:layers-1
    saliencyMaps2(c).mat = gaussProj2(c).Location - gaussProj2(c+1).Location;
    normMat2(c).mat = vecnorm((saliencyMaps2(c).mat)')';
end

%% Maximum Extraction 1

count=1;

figure
pcshow(ptCloud(1).Location,normMat(1).mat)
hold on

for c = 1:layers-1
for point = 1:size(ptCloud(1).Location,1)
    if ~ismember(point,wrongIndices1)
    [indices,dists] = findNeighborsInRadius(ptCloud(1),ptCloud(1).Location(point,:),((c))*baseSigma,'Sort',true);
    indices = setdiff(indices,wrongIndices1,'stable');
    flag = true;
    for i = 2:size(indices)
        if normMat(c).mat(indices(i)) >= normMat(c).mat(point)
            flag = false;
        end
    end
    if flag %&& size(indices,1) > 60*((c))
        maxIndices(count,c) = point;
        plot3(ptCloud(1).Location(point,1),ptCloud(1).Location(point,2),ptCloud(1).Location(point,3),'*r')
        count=count+1;
    end
    end
end
count=1;
end

plot3(ptCloud(1).Location(1,1),ptCloud(1).Location(1,2),ptCloud(1).Location(1,3),'*g')
hold off

%% Maximum extraction 2


count=1;

figure
pcshow(ptCloud(2).Location,normMat2(1).mat)
hold on

for c = 1:layers-1
for point = 1:size(ptCloud(2).Location,1)
    if ~ismember(point,wrongIndices2)
    [indices,dists] = findNeighborsInRadius(ptCloud(2),ptCloud(2).Location(point,:),((c))*baseSigma,'Sort',true);
    indices = setdiff(indices,wrongIndices2,'stable');
    flag = true;
    for i = 2:size(indices)
        if normMat2(c).mat(indices(i)) >= normMat2(c).mat(point)
           flag = false;
        end
    end
    if flag% && size(indices,1) > 60*(c)
        maxIndices2(count,c) = point;
        plot3(ptCloud(2).Location(point,1),ptCloud(2).Location(point,2),ptCloud(2).Location(point,3),'*r')
        count=count+1;
    end
    end
end
count=1;
end

plot3(ptCloud(2).Location(1,1),ptCloud(2).Location(1,2),ptCloud(2).Location(1,3),'*g')

hold off


%% Descriptor generation

M=4;
L=8;

saliencyNorms = zeros(size(normMat(1).mat,1),layers-1);

for c =1:layers-1
    for i = 1:size(normMat(1).mat,1)
        saliencyNorms(i,c)=normMat(c).mat(i);
    end
end

maxIndicesVector = reshape(maxIndices,[],1);
descriptorsVector=[]; 

missing = 0;

for c = 1:layers-1

for index=1:size(maxIndices,1)
    if maxIndices(index,c)>0
        descriptorsVector(index,c).zDir = normals(maxIndices(index,c),:);
        xyMat = null(descriptorsVector(index,c).zDir(:).');
        if abs(cross(xyMat(:,1),xyMat(:,2)) - descriptorsVector(index,c).zDir(:) ) < 0.0001
            descriptorsVector(index,c).xDir = xyMat(:,1);
            descriptorsVector(index,c).yDir = xyMat(:,2);
        else
            descriptorsVector(index,c).xDir = xyMat(:,2);
            descriptorsVector(index,c).yDir = xyMat(:,1);
        end

        [indicesOfNeigh,distances] = findNeighborsInRadius(ptCloud(1),ptCloud(1).Location(maxIndices(index,c),:),((c+1))*baseSigma,'Sort',true);
        neighborsGrid = zeros(size(indicesOfNeigh,1),3); %indices legend : M=1, L=2, CloudIndex=3
        for neighIndex = 1:size(indicesOfNeigh)
            diffVector = ptCloud(1).Location(maxIndices(index,c),:) - ptCloud(1).Location(indicesOfNeigh(neighIndex),:);      
            sigma = ((c+1))*baseSigma;
            neighborsGrid(neighIndex,1) = ceil((vecnorm(diffVector)*M)/sigma);
            phi = acos(dot(diffVector,descriptorsVector(index,c).zDir));
            normalOnPlane = (diffVector-(vecnorm(diffVector)*cos(phi) * descriptorsVector(index,c).zDir))/vecnorm(diffVector-(vecnorm(diffVector)*cos(phi) * descriptorsVector(index,c).zDir));
            normalsAndY=dot(normalOnPlane,descriptorsVector(index,c).yDir);
            if normalsAndY >= 0
                theta = acos(dot(normalOnPlane,descriptorsVector(index,c).xDir));
            else
                theta = 2*pi - acos(dot(normalOnPlane,descriptorsVector(index,c).xDir));
            end
            
            
            neighborsGrid(neighIndex,2) = ceil((theta*L)/(2*pi));
            neighborsGrid(neighIndex,3) = indicesOfNeigh(neighIndex);
        end
        neighborsGrid(neighborsGrid(:,1) > M,1) = M;
        neighborsGrid(neighborsGrid(:,1) < 1,1) = 1;
        neighborsGrid(isnan(neighborsGrid(:,2)),2) = 1;
        neighborsGrid = real(neighborsGrid);

        singleDescriptorGrid = zeros(M,L,2);
        singleDescriptorGridMask = ones(M,L,2);
        
        for range = 1:M
            for angle = 1:L
                %normalMean = mean(normals(M(:,1) == someValue, 3)); 
                  allM = neighborsGrid(:,1); % comma separated list expansion 
                  tfM = allM == range;
                  indexM = find(tfM) ;
                  allL = neighborsGrid(:,2); % comma separated list expansion 
                  tfL = allL == angle;
                  indexL = find(tfL) ;
                  
                  pointsInCell = intersect(indexM,indexL);
                  
                  averageNorm = [0,0,0];
                  averageSaliency = 0;
                  
                  if size(pointsInCell)>0
                    for i = 1:size(pointsInCell)
                        averageNorm = averageNorm + normals(neighborsGrid(pointsInCell(i),3),:);
                        averageSaliency = averageSaliency + saliencyNorms(neighborsGrid(pointsInCell(i),3),c);
                    end
                    averageNorm = averageNorm / norm(averageNorm);
                    averageSaliency = averageSaliency / size(pointsInCell,1);
                  else
                      singleDescriptorGridMask(range,angle,:) = [0 0];
                    
                  end
                  
                  
                  deltaNorm = 1.0 - abs(dot(averageNorm,normals(maxIndices(index,c),:)));
                  deltaSaliency = 1.0 - averageSaliency/saliencyNorms(maxIndices(index,c),c);
                  singleDescriptorGrid(range,angle,:) = [deltaNorm, deltaSaliency];
            end
        end
        
       % if maxIndicesVector(index)==23748 || maxIndicesVector(index)==25067 || maxIndicesVector(index)==5506
        maxIndicesVectorDescriptor(index-missing,c).descriptor = singleDescriptorGrid;
        maxIndicesVectorDescriptor(index-missing,c).mask = singleDescriptorGridMask;
        maxIndicesVectorDescriptor(index-missing,c).index = maxIndices(index,c);
       % else
      %      missing = missing+1;
      %  end
        
    end
end
end

%% Descriptor generation 2

saliencyNorms2 = zeros(size(saliencyMaps2(1).mat,1),layers-1);

for c =1:layers-1
    for i = 1:size(saliencyMaps2(1).mat,1)
        saliencyNorms2(i,c)=normMat2(c).mat(i);
    end
end

maxIndicesVector2 = reshape(maxIndices2,[],1);
descriptorsVector2=[]; 

missing = 0;
for c = 1:layers-1

for index=1:size(maxIndices2)
    if maxIndices2(index,c)>0
        descriptorsVector2(index,c).zDir = normals2(maxIndices2(index,c),:);
        xyMat = null(descriptorsVector2(index,c).zDir(:).');
        if abs(cross(xyMat(:,1),xyMat(:,2)) - descriptorsVector2(index,c).zDir(:) ) < 0.0001
            descriptorsVector2(index,c).xDir = xyMat(:,1);
            descriptorsVector2(index,c).yDir = xyMat(:,2);
        else
            descriptorsVector2(index,c).xDir = xyMat(:,2);
            descriptorsVector2(index,c).yDir = xyMat(:,1);
        end
        [indicesOfNeigh,distances] = findNeighborsInRadius(ptCloud(2),ptCloud(2).Location(maxIndices2(index,c),:),((c+1))*baseSigma,'Sort',true);
        neighborsGrid = zeros(size(indicesOfNeigh,1),3); %indices legend : M=1, L=2, CloudIndex=3
        for neighIndex = 1:size(indicesOfNeigh)
            diffVector = ptCloud(2).Location(maxIndices2(index,c),:) - ptCloud(2).Location(indicesOfNeigh(neighIndex),:);      
            sigma = ((c+1))*baseSigma;
            neighborsGrid(neighIndex,1) = ceil((vecnorm(diffVector)*M)/sigma);

            phi = acos(dot(diffVector,descriptorsVector2(index,c).zDir));
            normalOnPlane = (diffVector-(vecnorm(diffVector)*cos(phi) * descriptorsVector2(index,c).zDir))/vecnorm(diffVector-(vecnorm(diffVector)*cos(phi) * descriptorsVector2(index,c).zDir));
            normalsAndY=dot(normalOnPlane,descriptorsVector2(index,c).yDir);
            if normalsAndY >= 0
                theta = acos(dot(normalOnPlane,descriptorsVector2(index,c).xDir));
            else
                theta = 2*pi - acos(dot(normalOnPlane,descriptorsVector2(index,c).xDir));
            end
            
            neighborsGrid(neighIndex,2) = ceil((theta*L)/(2*pi));
            neighborsGrid(neighIndex,3) = indicesOfNeigh(neighIndex);
        end
        
        neighborsGrid(neighborsGrid(:,1) > M,1) = M;
        neighborsGrid(neighborsGrid(:,1) < 1,1) = 1;
        neighborsGrid(isnan(neighborsGrid(:,2)),2) = 1;
        neighborsGrid = real(neighborsGrid);
        
        singleDescriptorGrid2 = zeros(M,L,2);
        singleDescriptorGridMask2 = ones(M,L,2);

        
        for range = 1:M
            for angle = 1:L
                %normalMean = mean(normals(M(:,1) == someValue, 3)); 
                  allM = neighborsGrid(:,1); % comma separated list expansion 
                  tfM = allM == range;
                  indexM = find(tfM) ;
                  allL = neighborsGrid(:,2); % comma separated list expansion 
                  tfL = allL == angle;
                  indexL = find(tfL) ;
                  
                  pointsInCell = intersect(indexM,indexL);
                  
                  averageNorm = [0,0,0];
                  averageSaliency = 0;
                  
                  if size(pointsInCell)>0
                    for i = 1:size(pointsInCell)
                        averageNorm = averageNorm + normals2(neighborsGrid(pointsInCell(i),3),:);
                        averageSaliency = averageSaliency + saliencyNorms2(neighborsGrid(pointsInCell(i),3),c);
                    end
                    averageNorm = averageNorm / norm(averageNorm);
                    averageSaliency = averageSaliency / size(pointsInCell,1);
                  else
                      singleDescriptorGridMask2(range,angle,:) = [0 0];
                    
                  end
                  
                  deltaNorm = 1.0 - abs(dot(averageNorm,normals2(maxIndices2(index,c),:)));
                  deltaSaliency = 1.0 - averageSaliency/ saliencyNorms2(maxIndices2(index,c),c);
                  singleDescriptorGrid2(range,angle,:) = [deltaNorm, deltaSaliency];
            end
        end
        
       % if maxIndicesVector2(index)==10572 || maxIndicesVector2(index)==11891 || maxIndicesVector2(index)==2806
        maxIndicesVectorDescriptor2(index-missing,c).descriptor = singleDescriptorGrid2;
        maxIndicesVectorDescriptor2(index-missing,c).mask = singleDescriptorGridMask2;

        maxIndicesVectorDescriptor2(index-missing,c).index = maxIndices2(index,c);
       % else
       %     missing = missing +1;
       % end
        
    end
end
end


%% Feature Matching

length1=size(maxIndices,1);
length2=size(maxIndices2,1);
maxCompleteScore = zeros(size(maxIndicesVectorDescriptor,2),size(maxIndicesVectorDescriptor2,2));
maxCompleteScore = zeros(length1,length2,layers-1);

for c = 1:layers-1
for s = 1:length1
    for d = 1:length2
        maxCompleteScore(s,d,c) = 0;            
        if ~isempty(maxIndicesVectorDescriptor2(d,c).mask(:,:,1)) && ~isempty(maxIndicesVectorDescriptor(s,c).mask(:,:,1))
        for l_bar = 1:L
            normalScore = maxIndicesVectorDescriptor(s,c).mask(:,:,1).*circshift(maxIndicesVectorDescriptor2(d,c).mask(:,:,1),l_bar,2).*(1 - abs(maxIndicesVectorDescriptor(s,c).descriptor(:,:,1) - circshift(maxIndicesVectorDescriptor2(d,c).descriptor(:,:,1),l_bar,2)));
            saliencyScore = maxIndicesVectorDescriptor(s,c).mask(:,:,2).*circshift(maxIndicesVectorDescriptor2(d,c).mask(:,:,2),l_bar,2).*(1 - abs(maxIndicesVectorDescriptor(s,c).descriptor(:,:,2) - circshift(maxIndicesVectorDescriptor2(d,c).descriptor(:,:,2),l_bar,2)));
            
            normalScoreArray = normalScore';
            normalScoreArray = normalScoreArray(:)';
            
            saliencyScoreArray = saliencyScore';
            saliencyScoreArray = saliencyScoreArray(:)';
            
            tempCompleteSD = dot(normalScoreArray,saliencyScoreArray);
            if tempCompleteSD > maxCompleteScore(s,d,c)
                maxCompleteScore(s,d,c) = tempCompleteSD;
            end
        end
        end
        
    end
end
end

%% Selection of top 5 features


if iterator==9 || iterator==1
    numberOfSelected = 10;
else
    numberOfSelected = 5;
end

matchedFeatures = zeros(numberOfSelected,3);

for index=1:numberOfSelected
    [mxv,idx] = max(maxCompleteScore(:));
    [row,col,layer] = ind2sub(size(maxCompleteScore),idx);
    
    row = row(1);
    col = col(1);
    layer = layer(1);
    matchedFeatures(index,:) = [row col layer];
    maxCompleteScore(row,:,layer) = 0;
    maxCompleteScore(:,col,layer) = 0;
    if all(all(ismember(maxCompleteScore, max(maxCompleteScore(:))))) 
        break
    end

end


%% Get Point clouds from matched points

firstPoints = matchedFeatures(:,1) + (matchedFeatures(:,3)-1)*length1;
secondPoints = matchedFeatures(:,2) + (matchedFeatures(:,3)-1)*length2;


fixedCloud = ptCloud(1).Location(maxIndicesVector(firstPoints),:);
changingCloud = ptCloud(2).Location(maxIndicesVector2(secondPoints),:);



%% RANSAC

[tformEst,inlierIndex] = estimateGeometricTransform3D(changingCloud,fixedCloud,'rigid','maxNumTrials',100000,'Confidence',99); 

inliersPtCloud2 = pctransform(ptCloud(2),tformEst).Location;
inliersPtCloud1 = ptCloud(1).Location;

% inliersPtCloud2 = pctransform(pointCloud(changingCloud),tformEst).Location;
% inliersPtCloud1 = fixedCloud;


figure
firstPtCloud = pointCloud(inliersPtCloud1);
secondPtCloud = pointCloud(inliersPtCloud2);
pcshowpair(firstPtCloud,secondPtCloud)
title('Aligned point clouds')


%% ICP

ptCloudRef = firstPtCloud;
ptCloudCurrent = secondPtCloud;

gridSize = 0.001;
fixed = pcdownsample(ptCloudRef, 'gridAverage', gridSize);
moving = pcdownsample(ptCloudCurrent, 'gridAverage', gridSize);

tform = pcregistericp(moving, fixed, 'Metric','pointToPlane','Extrapolate', true);

ptCloudAlignedLocal = pctransform(ptCloudCurrent,tform);

ptCloudAligned = pctransform(ptCloudAlignedLocal,accumTform);
accumTform = rigid3d(accumTform.T * tformEst.T * tform.T );
mergeSize = 0.0015;
ptCloudScene = pcmerge(ptCloudScene, ptCloudAligned, mergeSize);
% Visualize the world scene.

ptLocal = pcmerge(ptCloudRef, ptCloudAlignedLocal, mergeSize);


figure
pcshow(ptLocal, 'VerticalAxis','Y', 'VerticalAxisDir', 'Down')
drawnow



figure
hAxes = pcshow(ptCloudScene, 'VerticalAxis','Y', 'VerticalAxisDir', 'Down');
title('Updated world scene')
% Set the axes property for faster rendering
hAxes.CameraViewAngleMode = 'auto';
hScatter = hAxes.Children;
drawnow    
end

function smoothedPoint = smoothPointInCloud(ptCloud, point,sigma)
    [neighbors, dists] = findNeighborsInRadius(ptCloud,point,2*sigma,'Sort',true);
    neighbors = ptCloud.Location(neighbors,:);
    denom = sum( exp( -((dists.^2)/(2*sigma^2) )) );
    numer = sum(neighbors.* exp( -((dists.^2)/(2*sigma^2) )) ,1);
    smoothedPoint = numer/denom;
end

function centeredCloud = centerCloud(ptCloud)
    xyz0(1)=(ptCloud.XLimits(1)+ptCloud.XLimits(2))/2;
    xyz0(2)=(ptCloud.YLimits(1)+ptCloud.YLimits(2))/2;
    xyz0(3)=(ptCloud.ZLimits(1)+ptCloud.ZLimits(2))/2;

    centeredCloud=pointCloud(bsxfun(@minus,ptCloud.Location,xyz0));
end

function wrongIndices= cleanCloud(ptCloud,baseSigma)
    locations = [];
    wrongIndices = [];
    counter = 1;
    counter2 = 1;
    for point = 1:ptCloud.Count
        [indices,dists] = findNeighborsInRadius(ptCloud,ptCloud.Location(point,:),baseSigma,'Sort',true);
        neighborhood = pointCloud(ptCloud.Location(indices,:));
        midPoint = ([neighborhood.XLimits(1) neighborhood.YLimits(1) neighborhood.ZLimits(1)] + [neighborhood.XLimits(2) neighborhood.YLimits(2) neighborhood.ZLimits(2)])./2;
        span = norm([neighborhood.XLimits(1) neighborhood.YLimits(1) neighborhood.ZLimits(1)] - [neighborhood.XLimits(2) neighborhood.YLimits(2) neighborhood.ZLimits(2)]);
        if norm(ptCloud.Location(point,:) - midPoint)<=span/6
            locations(counter) = point;
            counter = counter+1;
        else
            wrongIndices(counter2) = point;
            counter2 = counter2+1;

        end
    end
    cleanedCloud=pointCloud(ptCloud.Location(locations,:));
    
end

