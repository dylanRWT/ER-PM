%% LAMAprep
% Chooses a folder, detects all .txt files and performs conversions on them
% so only frame,x,y,intensity are imported for later use in LAMA
%
inputdir = uigetdir();  % choose the directory where the .csv files are held
files=dir(strcat(inputdir,'/*.csv')); 
files={files.name};
outputdir = uigetdir(inputdir);     % choose where you want the new.txt files for lama to go
%%
dirlist = {};
for file = 1:width(files)
    temptable = readtable(strcat(inputdir,'/',files{file}));
    newmatrix = [];
    newmatrix = [temptable.x_nm_ , temptable.y_nm_, temptable.frame, temptable.intensity_photon_];
    writematrix(newmatrix,strcat(outputdir,'/',extractBefore(files{file},'_full.csv')))
    dirlist{file,1} = strcat(outputdir,'/',extractBefore(files{file},'_full.csv'),'.txt');
end
writecell(dirlist,strcat(outputdir,'/registry.txt')) % make a registry so you can quickly process the LAMA files in LAMA
