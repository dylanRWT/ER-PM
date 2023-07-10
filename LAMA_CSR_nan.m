%% LAMA_CSR_nan
%
% Created by Dylan Greening for use in RWTurner lab UCalgary 
%
% Input Arg:
% indicies = all indicies which need to be pooled in the final table
% array1 = an array where each cell is a matrix of nnd values 
% 
% Output Arg:
% matrixout = a matrix where each column is a diffrent data set and height
%   differences between columns are filled with nan
%
%%
function matrixout = LAMA_CSR_nan(array1,indicies)
% a function which takes arrays of cell nnd and reorganizes them into a
% single table using nan to fill in the gaps.
% this is done for later assessment by gammMixtureMaker_singleimages 
%
for i = 1:height(indicies) %look at ever index pair for the cells
    if i == 1
    matrixout(:,i) = array1{indicies(i,1),indicies(i,2)}; 
    maxheight = height(array1{indicies(i,1),indicies(i,2)});
    else 
        tempheight = height(array1{indicies(i,1),indicies(i,2)});
        if tempheight < maxheight
            diffheight = maxheight-tempheight;
            tempdata = [array1{indicies(i,1),indicies(i,2)};nan(diffheight,1)];
            matrixout(:,i) = tempdata;
        elseif tempheight > maxheight
            diffheight = tempheight - maxheight;
            maxheight = tempheight;
            matrixout = [matrixout; nan(diffheight,width(matrixout))];
            matrixout(:,i) = array1{indicies(i,1),indicies(i,2)};
        elseif tempheight == maxheight
            matrixout(:,i) = array1{indicies(i,1),indicies(i,2)};
        end
    end
end
