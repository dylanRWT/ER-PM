%% UIA_NNMF
% function used by runRefineOnOldWorkspaces.m
% Written by Dylan Greening for RWTurner Lab UCalgary
function [nnmf_output,centroidCounts,centroidPositions] = UIA_NNMF(img_array,mask_array,thresh,px,py,numOfFeatures,SegPxRad)
    nnmf_output = []; % empty array to become .csv
    centroidCounts = []; % empty array to hold centroid counts
    filenum = width(img_array); 
    centroidPositions = {};
    for i = 1:filenum % for each file in the workspace
       strcat("processing ", num2str(i), " of ", num2str(filenum)) % this counts down so you get less impatient with the process
        centroidX = UIA_centroid_refined((img_array{i}.*mask_array{i}),thresh,SegPxRad); % find all the centroids
        centroidPositions{i} = centroidX; 
        [~,~,rotrawX] = UIA_img_rotate(img_array{i}.*mask_array{i},centroidX,thresh,px,py); % generate centroid images with Nearest Neighbor to the left
        [~,h] = nnmf(rotrawX,numOfFeatures);    % plug all these centroid images into non-negative matrix factorization
        centroidCounts(end+1) = height(centroidX);          % store the number of centroids detected in this one image
        for j = 1:numOfFeatures                 % reshape the linear centroid images into meaningful matricies
            temp = reshape(h(j,:),px*2,py*2);   
            A = temp(px,:);                     % take the middle row 
            nnmf_output(:,end+1) = A';          % store it in nnmf_output
        end
    end  
end
