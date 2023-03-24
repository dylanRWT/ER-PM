% UIA - Unsupervised Image Analysis
% A script concieved by Wilten Nicola for use in the RW Turner Lab
% Altered for specific applications by Dylan Greening 
% Takes one or two channel *.tif files masks the soma and generates an
% overlap image 
% For use at the University of Calgary.
% Last modifed 2023 march 01 
% Requiers functions : 
% UIA_imloader_ssV2 
% UIA_imoverlap_v3 
%
% Input 
% Folder containing files ending in "_568.tif" or "_647.tif" and an output
% folder where .mat files will be saved
clear all
close all
clc
%% make a list of all files in the directory
inputdir = "*/Input folder/" ;
outputdir = "*/Output folder/" ;
%%
mixeddata = 1;
if mixeddata == 0    
    file_struct_568 = dir(strcat(inputdir,"*_568.tif"));    % get a structure of the directory
    file_list1 = {file_struct_568(:).name} ;                % get just file names as a cell array
    file_struct_647 = dir(strcat(inputdir,'*_647.tif'));    % get a structure of the directory
    file_list2 = {file_struct_647(:).name} ;                % get just file names as a cell array
else 
    file_struct_X = dir(strcat(inputdir,'*.tif'));          % get a structure of the directory
    file_listX = {file_struct_X(:).name} ;                  % get just file names as a cell array
end
px = 300;   % this is half of the size of the centroid images for later
py = 300;
% initialize data frames and arrays
if mixeddata == 0 
    filenum = width(file_list2);
    img_array_red = {};
    img_array_green = {};
elseif mixeddata == 1
    filenum = width(file_listX);
    img_array = {};
end
mask_array = {};
% open each image and store its masking positions
for i = 1:filenum
    if mixeddata == 0 
        currfile2 = strcat(inputdir,file_list2{i});
        currfile1 = strcat(inputdir,file_list1{i});
        [img_array_red{i},mask_array{i}] = UIA_imloader_ssv2(currfile2,1,i,filenum);     
        [img_array_green{i},~] = UIA_imloader_ssv2(currfile1,0,i,filenum); 
    elseif mixeddata == 1
        currfile = strcat(inputdir,file_listX{i});
        [img_array{i},mask_array{i}] = UIA_imloader_ssv2(currfile,1,i,filenum);
    else
        error("there was no file list provided for processing as an image")
    end
end
prompt = " What name do you want to use for this dataset? ";
lil_nam = inputdlg(prompt);
% Generate overlap images; 
close all
if mixeddata == 0 && (height(file_struct_647)>0 && height(file_struct_568)>0)
    img_overlap = {};
    for i =1:filenum
        [img_overlap{i}] = UIA_imoverlap_v3({img_array_red{i},img_array_green{i}});
        temp = extractBefore(file_list1{i},"_568");
        imwrite(img_overlap{i},strcat(outputdir,temp,"_overlap.tif"))
    end
end
save(strcat(outputdir,lil_nam,"gauss_workspace.mat"))