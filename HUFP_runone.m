%% HUPF_runone  
%    
% Created by Dylan Greening for RW Turner Lab UCalgary
% with input from Wilten Nicola of UCalgary
%
% Created March 28 2022
% Last modified 2023 Jan
% 
% Requiers Scripts:
% HUFP_peak_detection
% HUFP_ladder_analysis
% HUFP_blob_analysis
% HUFP_residual_plotter
% HUFP_osscilation_fitter_auto
% HUFP_osscilation_fitter_manual
% 
% Selects a .csv of interest produced from NNMF of RGB-tiffs
% finds the displacement of peaks and plots said peaks to identify
% periodicity in biological structures. 
clear all
clc
% PARAMETERS
folder_name = '../HUFP output/';
feat_per_cell = 6; % number of features per image
px_sz = 7.52     ; % true size of one subpixel
x0 = 1:600      ;  % initialized x-variable
% Open .csv
[file,path] = uigetfile('../HUFP input/*.csv') ; % select the file
    x_channel = 0;
    lwcp = 15         ; % lower crop range for x axis data of peak distances usually 15
    upcp = lwcp + 55   ; % upper crop range for x axis data of peak distances usually 55
raw = readmatrix(strcat(path,file)) ;     % read the .csv
base_file_name = extractBefore(file,'_') ;  % make a variable to use save info in 
feat_num = width(raw); % 
% Build a table for peak heights to be entered into
[sm_fm,pk_fm,ladder_feats,blob_feats,neither_feats,blob_count] = HUFP_peak_detection(x0,raw,feat_num,feat_per_cell) ;
% Evaluate ladder classified features
jam = input("Do you want to evaluate ladders? y or n  ", 's') ;
close all
if jam == 'y'
    [best_period_px,decent_displacements,~,~,ivt,best_model,displacements] = HUFP_ladder_analysis(x0,pk_fm,ladder_feats,lwcp,upcp, folder_name, base_file_name, px_sz) ;
end
% Evaluate blob classified features positions
jam = input("Do you want to evaluate blobs? y or n ", 's') ;
if jam == 'y'
    HUFP_blob_analysis(blob_count,blob_feats,sm_fm,folder_name,base_file_name,px_sz)
end
%
save(strcat(folder_name,"matlab workspaces/",base_file_name,".mat"))
%
disp(round(best_period_px*px_sz,8))