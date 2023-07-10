%% LAMA_CSR1
%
% Created by Dylan Greening for use in RWTurner Lab UCALGARY
%
% A function used in LAMA_CSR.m
%
% Uses centroid positions and binary image masks to define which centroids
% fall without the mask
%
% Input Arg:
% mca_positions = array where each cell defines a nx2 matrix referenceing
%   (x,y) coordinates of centroids. Cells in the first column are for 647
%   channel while the second column is the 568 channel. 
% binaryImages = an array where each cell is a square binary image which
%   acts as a mask
%
% Output Arg:
% mca_culled647 = array with cells containing a nx2 matrix where rows are
%   (x,y) coordinates of in mask centroids. 647 is assumed the first channel
%   in this workflow
% mca_culled568 = like mca_culled647 but for the second channel
% sum_of_centroidsXXX = number of centroids detected in the mask
% tango = unk
%
%
function [mca_culled647,sum_of_centroids647,mca_culled568,sum_of_centroids568,tango] = LAMA_CSR1(mca_positions,binaryImages)
    %%
    mca_culled647 = {};
    sum_of_centroids647 = [];
    mca_culled568 = {};
    sum_of_centroids568 = [];
    tango = {};
    %
    mca_position_in_mask = {};
    for i = 1:height(mca_positions) % look at each cell in the first (647) column and assess which centroids fall within the mask
        if find(~mca_positions{i,1}(:,5))          % This is a small correction for the edge case where the rounding of a centroid position results in a 0 px value
            I = find(~mca_positions{i,1}(:,5));                % and corrects this by providing a value of 1
            for w = 1:height(I)                         
               mca_positions{i,1}(I(w),5) = 1;
            end
        end
        if find(~mca_positions{i,1}(:,6))
            I = find(~mca_positions{i,1}(:,6));
            for w = 1:height(I)
                mca_positions{i,1}(I(w),6) = 1;
            end
        end
        mca_position_in_mask{i,1} = zeros(height(mca_positions),1); % initialize a matrix of zeros the same height as the centroids
        for j = 1:height(mca_positions{i})                                           
            if binaryImages{i}(mca_positions{i,1}(j,6),mca_positions{i,1}(j,5)) % check if the pixel position is within the mask
                mca_position_in_mask{i,1}(j,1) = 1;
            end
        end
        sum_of_centroids647(i) = sum(mca_position_in_mask{i}); % store the number of centroids in the mask for randomization later
        mca_culled647{i,1} = mca_positions{i,1}(logical(mca_position_in_mask{i,1}),:); % use the logical index to cull the original positions to only thoes in the mask
        tango{i,1} = all(ismember(mca_culled647{i,1},mca_positions{i,1}),2);
    end
    mca_position_in_mask = {};
    for i = 1:height(mca_positions)
            if find(~mca_positions{i,2}(:,5))
            I = find(~mca_positions{i,2}(:,5));
            for w = 1:height(I)
                mca_positions{i,2}(I(w),5) = 1;
            end
        end
        if find(~mca_positions{i,2}(:,6))
            I = find(~mca_positions{i,2}(:,6));
            for w = 1:height(I)
                mca_positions{i,2}(I(w),6) = 1;
            end
        end
        mca_position_in_mask{i,1} = zeros(height(mca_positions),1); 
        for j = 1:height(mca_positions{i,2})
            if binaryImages{i,1}(mca_positions{i,2}(j,6),mca_positions{i,2}(j,5)) % check if the pixel position is within the mask
                mca_position_in_mask{i,1}(j,1) = 1;
            end
        end
        sum_of_centroids568(i) = sum(mca_position_in_mask{i,1}); % store the number of centroids in the mask for randomization later
        mca_culled568{i,1} = mca_positions{i,2}(logical(mca_position_in_mask{i,1}),:); % use the logical index to cull the original positions to only thoes in the mask
        tango{i,2} = all(ismember(mca_culled568{i,1},mca_positions{i,2}),2);
    end
    %%
end