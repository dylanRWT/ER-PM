% UIA_centroid_refined
% function for use with runRefineOnOldWorkspaces.m and nested within
% UIA_NNMF.m
%
% function performes a segmentation of the image based off of grayscale intensity
%  and returns the weighted centroids for the image
function centroids = UIA_centroid_refined(img,thresh,SegmentationPxRadius)
    se = strel('disk',SegmentationPxRadius);   % this segments the image into layers of the prescribed shape
    img_segmented = imopen(img,se); % this opens the image using the segmentation information
    %immax = imregionalmax(img_segmented); % this works to find regional maxima 
    weightedCentroidPositions = regionprops(img.*imregionalmax(img_segmented)>thresh,img,'WeightedCentroid'); 
    centroids = zeros(height(weightedCentroidPositions),2);
    for k = 1:height(weightedCentroidPositions)
    centroids(k,:) = [weightedCentroidPositions(k).WeightedCentroid(1),weightedCentroidPositions(k).WeightedCentroid(2)];
    end
end