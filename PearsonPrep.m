%% PearsonPrep
% A script for converting images from RBG tiff into 8-bit gray scale and
% masked for soma only from binary images
imgdir = strcat(uigetdir('Select a folder with all the images in it '),'/');
maskdir = strcat(uigetdir('Select a folder with the relevant binary images '),'/');
outputdir = strcat(uigetdir('Select an output folder '),'/');
%%
imglist = dir(strcat(imgdir,"*.tif"));
imglist = {imglist.name};
masklist = dir(strcat(maskdir,"*.tif"));
masklist = {masklist.name};
% look at each image
for i = 1:width(imglist)
    currentimage = rgb2gray(imread(strcat(imgdir,imglist{i})));
    tempname = extractBefore(imglist{i},'_');
    for j = 1:width(masklist)
        if strcmp(tempname,extractBefore(masklist{j},'_'))
            tempmask = imread(strcat(maskdir,masklist{j}));
        end
    end
    maskedimage = currentimage .* uint8(tempmask);
    imwrite(maskedimage,strcat(outputdir,extractBefore(imglist{i},'.tif'),'mask.tif'))
end
