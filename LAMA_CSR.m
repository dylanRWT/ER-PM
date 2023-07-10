%% LAMA_CSR
%
% Created by Dylan Greening for use in RWTurner Lab UCALGARY
%
% CSR of LAMA defined Centroids
% Takes mca_analysis.txt files generated in LAMA masks the images with a
% somatic binary image defined by the NNMF workflow, defined
% nearest-neighbour distances in that euclidian space and recalculates with
% randomized positions.
% 
% Single Symbol Key: First symbol refers to 647 channel second symbol the 568 channel
% 1 = actinin1
% 2 = actinin2
% C = Cav1.3
% R = RyR2
% K = IK (KCa3.1)
% S = SpectrinBii
% N = Nav1.6
% B = AnkyrinB
%
% Input:
% mca_analysis.txt
% pixel magnification
% binary image
%
% Output:
% NND.csv
% ANN 
% csv of NND (for input through gammaMixtureMaker)
% ripleysK 
% radial distribution function plots
%
% Functions:
% LAMA_CSR1 - initializes mask information
% LAMA_CSR2 - performs the randomization of the centroids in the mask
% LAMA_CSR3 - performs the statistics on the centroid positions
% LAMA_CSR_ann - script to calc average nearest neighbour
% LAMA_CSR_distrange - script to calc radial disribution of datasets
% LAMA_CSR_err - script to calc mean stat and SEM of mean
% LAMA_CSR_nan - script to format NND into a matrix
% LAMA_CSR_nnd1 - script to pull nearest neighbour distributions  
% LAMA_CSR_pdist1 - small script to perform repetitive pdist calcs
%
clear all
clc
close all
%
debugger = 0;
visualize = 0;
num_of_runs = 1; % the number of position randomizations to perform per image
px_mag = 20; % subpixelization in STORM reconstruction
px_size = 150.5; % pixel size in nanometers on raw camera
nm_to_px = px_mag/px_size;
master_dir = uigetdir(".."); % select the parent folder with mca_analysis files in it
lil_nam = master_dir(end-1:end);
%
mca_dir_files = dir(strcat(master_dir,"/",lil_nam,"mca/*647.txt")); % pull only .txt files
mca_dir_files = {mca_dir_files.name}; % and then only store their names
mca_positions = {}; % a new variable to only hold the centroid positions
%
for i = 1:width(mca_dir_files)
    temp = readmatrix(strcat(master_dir,"/",lil_nam,"mca/",mca_dir_files{i}));
    mca_positions{i,1} = [temp(:,1:2),temp(:,1:2)*nm_to_px,round(temp(:,1:2)*nm_to_px)];   
end
%
mca_dir_files = dir(strcat(master_dir,"/",lil_nam,"mca/*568.txt")); % pull other channel files
mca_dir_files = {mca_dir_files.name}; % and then only store their names
%
for i = 1:width(mca_dir_files)
    temp = readmatrix(strcat(master_dir,"/",lil_nam,"mca/",mca_dir_files{i}));
    mca_positions{i,2} = [temp(:,1:2),temp(:,1:2)*nm_to_px,round(temp(:,1:2)*nm_to_px)];     
end
%
bi_dir_files = dir(strcat(master_dir,"/",lil_nam,"_binarymasks/*mask.tif")); % only collect names of the masks
bi_dir_files = {bi_dir_files.name};    
binaryImages = {};      %initialize a variable to hold the masks
%
for i = 1:width(bi_dir_files)
    binaryImages{i,1} = imread(strcat(master_dir,"/",lil_nam,"_binarymasks/",bi_dir_files{i}));
end
% convert the mca_positions to pixel distances temporarily
bi_im_limits = zeros(1,4);
for i = 1:height(mca_positions)
    minX = width(binaryImages{i,1}) + 1;
    maxX = 0;
    minY = height(binaryImages{i,1})+ 1;
    maxY = 0;
    for row = 1:width(binaryImages{i})
        if sum(find(binaryImages{i,1}(row,:))) > 0 
            if find(binaryImages{i}(row,:),1,'first') < minX 
                minX = find(binaryImages{i}(row,:),1,'first');
            end
            if find(binaryImages{i}(row,:),1,'last') > maxX
                maxX = find(binaryImages{i}(row,:),1,'last');
            end
        end
    end
    for col = 1:height(binaryImages{i})
        if sum(find(binaryImages{i}(:,col))) > 0 
            if find(binaryImages{i}(:,col),1,'first') < minY
                minY = find(binaryImages{i}(:,col),1,'first');
            end
            if find(binaryImages{i}(:,col),1,'last') > maxY
                maxY = find(binaryImages{i}(:,col),1,'last');
            end
        end
    end
    bi_im_limits(i,:) = [minX,maxX,minY,maxY]; % minX maxX minY maxY
end
%
if debugger == 1
    for i = 1:height(mca_positions)
        figure()
        imshow(binaryImages{i}), hold on
        scatter(mca_positions{i,1}(:,3),mca_positions{i,1}(:,4),'.g')
        scatter(mca_positions{i,2}(:,3),mca_positions{i,2}(:,4),'.b')
        xline(bi_im_limits(i,1),'r')
        xline(bi_im_limits(i,2),'r')
        yline(bi_im_limits(i,3),'r')
        yline(bi_im_limits(i,4),'r'), hold off
    end
    %
    uinput = input("Enter 1 if the mca & mask files seem aligned     ");
    if uinput == 1
        close all
    else  
        error(" You mentioned that the mca_analysis and binary masks were not aligned, go fix this ")
    end
end
% Find masked centroids
% return positions and counts and logic position for what was counted in
% from mca_positions 
[mca_culled647,sum_of_centroids647,mca_culled568,sum_of_centroids568,tango] = LAMA_CSR1(mca_positions,binaryImages); % 
% Randomize the centroids in the masks
csr_cent_647 = LAMA_CSR2(binaryImages,bi_im_limits,sum_of_centroids647,num_of_runs,nm_to_px); % for the 647 channel (first coded letter)
csr_cent_568 = LAMA_CSR2(binaryImages,bi_im_limits,sum_of_centroids568,num_of_runs,nm_to_px); % and the 568 channel (second coded letter)
if debugger ==1 
    % show the masks with the random positions to ensure everything went ok
    for i = 1:height(mca_positions) % only show on column
        figure()
        imshow(binaryImages{i}), hold on
        scatter(csr_cent_647{i,1}(:,1).*nm_to_px,csr_cent_647{i,1}(:,2).*nm_to_px,'.g')
        scatter(csr_cent_568{i,1}(:,1).*nm_to_px,csr_cent_568{i,1}(:,2).*nm_to_px,'.b')
        xline(bi_im_limits(i,1),'r')
        xline(bi_im_limits(i,2),'r')
        yline(bi_im_limits(i,3),'r')
        yline(bi_im_limits(i,4),'r'), hold off
    end
    uinput = input("Enter 1 if the mca & mask files seem aligned     ");
    if uinput == 1
        close all
    else  
        error(" You mentioned that the mca_analysis and binary masks were not aligned, go fix this ")
    end
end
%% run stats
[nnd_singles,h] = LAMA_CSR3(mca_culled647,csr_cent_647,mca_culled568,csr_cent_568,0:10:2000,visualize,lil_nam);
%
mkdir(strcat(master_dir,'/NND/'))
if visualize == 1
    savefig(h,strcat(master_dir,'/NND/',lil_nam,'figures.fig'))
end
%
close all
writematrix(nnd_singles{1,1},strcat(master_dir,'/NND/',lil_nam(1),lil_nam(1),'_NND','.csv'))
writematrix(nnd_singles{1,2},strcat(master_dir,'/NND/',lil_nam(2),lil_nam(2),'_NND','.csv'))
writematrix(nnd_singles{1,3},strcat(master_dir,'/NND/',lil_nam(1),lil_nam(1),'csr_NND','.csv'))
writematrix(nnd_singles{1,4},strcat(master_dir,'/NND/',lil_nam(2),lil_nam(2),'csr_NND','.csv'))
writematrix(nnd_singles{1,5},strcat(master_dir,'/NND/',lil_nam(1),lil_nam(2),'_NND','.csv'))
writematrix(nnd_singles{1,6},strcat(master_dir,'/NND/',lil_nam(1),lil_nam(2),'csr_NND','.csv'))
%



