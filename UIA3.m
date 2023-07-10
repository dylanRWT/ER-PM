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
inputdir = strcat(uigetdir('Select an input folder'),'/');
outputdir = strcat(uigetdir('Select an output folder'),'/');
%% 
file_struct_568 = dir(strcat(inputdir,"*_568.tif"));    % get a structure of the directory
file_list1 = {file_struct_568(:).name} ;                % get just file names as a cell array
file_struct_647 = dir(strcat(inputdir,'*_647.tif'));    % get a structure of the directory
file_list2 = {file_struct_647(:).name} ;                % get just file names as a cell array
px = 300;   % this is half of the size of the centroid images for later
py = 300;
% initialize data frames and arrays
filenum = width(file_list2);
img_array_red = {};
img_array_green = {};
mask_array = {};
% open each image and store its masking positions
for i = 1:filenum
        currfile2 = strcat(inputdir,file_list2{i});
        currfile1 = strcat(inputdir,file_list1{i});
        [img_array_red{i},mask_array{i}] = UIA_imloader_ssv2(currfile2,1,i,filenum);     
        [img_array_green{i},~] = UIA_imloader_ssv2(currfile1,0,i,filenum); 
end
prompt = " What name do you want to use for this dataset? ";
lil_nam = inputdlg(prompt);
%% Generate overlap & mask images;
%
outputdir = strcat(outputdir,"/",lil_nam);
mkdir(outputdir);
close all
outputA = strcat(outputdir,"/",lil_nam,"_overlaps");
outputB = strcat(outputdir,"/",lil_nam,"_binarymasks");
mkdir(outputA)
mkdir(outputB)
img_overlap = {};
for i =1:filenum
    [img_overlap{i}] = UIA_imoverlap_v3({img_array_red{i},img_array_green{i}});
    temp = extractBefore(file_list1{i},"_568");
    imwrite(img_overlap{i},strcat(outputA,"/",temp,"_overlap.tif"))
    imwrite(mask_array{i},strcat(outputB,"/",temp,"_mask.tif"))
end
%% Conduct NNMF on the masked images
thresh = 0.245; % this is a threshold for centroid detection
diskSize = 2;   % the size of disks to use when segmenting the images for the weighted centroid detection
nof = 6;        % this is the number of features pulled out of each masked parent image
[XRed,hcRed,centRed] = UIA_NNMF(img_array_red,mask_array,thresh,px,px,nof,diskSize); 
writematrix(XRed,strcat(outputdir,"/",lil_nam,"_refinedR.csv")) 
[XGreen,hcGreen,centGreen] = UIA_NNMF(img_array_green,mask_array,thresh,px,px,nof,diskSize); 
writematrix(XGreen,strcat(outputdir,"/",lil_nam,"_refinedG.csv")) 
[XOverlap,hcOverlap,centOverlap] = UIA_NNMF(img_overlap,mask_array,thresh,px,px,nof,diskSize); 
writematrix(XOverlap,strcat(outputdir,"/",lil_nam,"_refinedOverlap.csv")) 
writematrix(XOverlap,strcat("../HUFP input/",lil_nam,num2str(nof),"_refinedOverlap.csv"))
%
save(strcat(outputdir,"/",lil_nam,"_gauss_workspace.mat"))
%%
figure()
ax(1) = subplot(1,2,1) ;
imshow(img_array_green{4})
ax(2) = subplot(1,2,2);
imshow(img_array_green{4}),hold on
scatter(centGreen{4}(:,1),centGreen{4}(:,2),'xr'), hold off
linkaxes(ax)