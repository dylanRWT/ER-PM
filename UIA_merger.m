%% UIA_merger
% A function which uses workspaces genereated by the UIA3 merging the data
% from each individual workspace into a larger dataset for use in
% HUFP_runone
clear all
clc
close all
%
dirnames = dir("../UIA output");
dirnames = struct2table(dirnames);
dirnames = dirnames.name(dirnames.isdir);
dirnames = dirnames(~ismember(dirnames,{'.','..'}));
outputfolder = "../HUFP input";
% initialize possible feature sets
NN = zeros(600,1);
a11 = NN;
a22 = NN;
SS = NN;
BB = NN;
KK = NN;
CC = NN;
RR = NN;
% Look at each folder and store the appropriate data into the initialized feature sets
for folder = 1:height(dirnames)
    smallfolder = dir(strcat("../UIA output/",dirnames{folder,1},"/*.csv"));
    smallfolder = {smallfolder.name};
    for file = 1:width(smallfolder)
        tempmatrix = [];
        temp = extractBefore(smallfolder{1,file},".csv");
        if strcmp(temp(end),"G")
            tempnum = 2;
            tempmatrix = readmatrix(strcat("../UIA output/",dirnames{folder,1},'/',smallfolder{1,file}));
        elseif strcmp(temp(end),"R")
            tempnum = 1;
            tempmatrix = readmatrix(strcat("../UIA output/",dirnames{folder,1},'/',smallfolder{1,file}));
        end
        if all(size(tempmatrix) ~= [0,0]);
            if  strcmp(temp(tempnum),'N')
                NN = [NN,tempmatrix];
            elseif  strcmp(temp(tempnum),'C')
                CC = [CC,tempmatrix];
            elseif  strcmp(temp(tempnum),'B')
                BB = [BB,tempmatrix];
            elseif  strcmp(temp(tempnum),'K')
                KK = [KK,tempmatrix];
            elseif  strcmp(temp(tempnum),'R')
                RR = [RR,tempmatrix];
            elseif  strcmp(temp(tempnum),'S')
                SS = [SS,tempmatrix];
            elseif  strcmp(temp(tempnum),'1')
                a11 = [a11,tempmatrix];
            elseif  strcmp(temp(tempnum),'2')
                a22 = [a22,tempmatrix];
            end
        end
    end
end
% remove unnecessary column used only for initialization
rmcol = sum(NN,1) == 0;
NN = NN(:,~rmcol);
%
rmcol = sum(CC,1) == 0;
CC = CC(:,~rmcol);
%
rmcol = sum(KK,1) == 0;
KK = KK(:,~rmcol);
%
rmcol = sum(BB,1) == 0;
BB = BB(:,~rmcol);
%
rmcol = sum(RR,1) == 0;
RR = RR(:,~rmcol);
%
rmcol = sum(SS,1) == 0;
SS = SS(:,~rmcol);
%
rmcol = sum(a11,1) == 0;
a11 = a11(:,~rmcol);
%
rmcol = sum(a22,1) == 0;
a22 = a22(:,~rmcol);
% save .csv files
writematrix(NN,strcat(outputfolder,'/NN6_refined.csv'))
writematrix(BB,strcat(outputfolder,'/BB6_refined.csv'))
writematrix(CC,strcat(outputfolder,'/CC6_refined.csv'))
writematrix(KK,strcat(outputfolder,'/KK6_refined.csv'))
writematrix(RR,strcat(outputfolder,'/RR6_refined.csv'))
writematrix(SS,strcat(outputfolder,'/SS6_refined.csv'))
writematrix(a11,strcat(outputfolder,'/116_refined.csv'))
writematrix(a22,strcat(outputfolder,'/226_refined.csv'))

