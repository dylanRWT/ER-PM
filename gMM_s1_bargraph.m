%% gMM_s1_bargraph
%
% Created by Dylan Greening for use in RWTurner Lab UCALGARY
%
% function for use in gammaMixtureMaker_singles
% 
% takes data from the mastertable generated in gammaMixtureMaker_singles
% and produces bar graphs with SEM for visual comparison.
%
% Input arg:
% mastertable = a table of row > 1 and col = 22
% setNames = an array of strings which represent groupings in the
%   mastertable
% CODlim = R2 value which acts as a cutoff for data incorporation
% mixlim = mixing portion cutoff when peaks == 2
% peaks = single gamma fit or double gamma fit
% visualize = logical input for graphics generation
% stat = a string either 'median','mean', or 'mode'
%
% Output arguments:
% mu = a 2 row, x col array where the upper row represents the mean stats
%   of the near zero gamma distribution and the lower row represents the mean
%   stats of the far from zero gamma distribution
% sem = like mu but representing the standard error of the mean of the near
%   or far gamma mixture model stats
%%
function [mu,sem,output] = gMM_s1_bargraph(mastertable,setNames,CODlim,mixlim,peaks,visualize,stat)
% determine which stat to visualize and report
if strcmp(stat,'median')
    stat1 = [7,15];
elseif strcmp(stat,'mean')
    stat1 = [8,16];
elseif strcmp(stat,'mode')
    stat1 = [6,14];
end
output={};
nstat = {}; % near stats
fstat = {}; % far stats
uptemp = table2cell(mastertable); % convert the mastertable to a easy to work with cell structure
if peaks == 1 % if there is only one fit distribution as defined by the user
    for set = 1:width(setNames) % look up every grouping
       uselist = false(height(uptemp),1); % and generate a logical array to define useful entries for that grouping
        for i = 1:height(uptemp)            
            if strcmp(uptemp{i,2},setNames{set})
                if uptemp{i,5} > CODlim  % cull based off of COD                    
                        uselist(i,1) = 1;
                end
            end
        end
        tempstat{1,set} = cell2mat(uptemp(uselist,stat1(1))); % if the entry survived filtering add it to a temparray
        output{1,set} = tempstat{1,set};
    end
   
    if visualize == 1
        figure()
            [mu,~,~,sem] = LAMA_CSR_err(tempstat); % use this function to find the error in the 
            bar(1:width(setNames),mu) , hold on 
            er = errorbar(1:width(setNames),mu,sem,sem);
            er.Color = [0 0 0];                            
            er.LineStyle = 'none'; 
            hold off
            ylim([0 1000])
            xticklabels(extractBefore(setNames,'_NND'))
            ylabel(strcat("Near ",stat," (nm)"))
    else
        [mu,~,~,sem] = LAMA_CSR_err(tempstat);
    end
elseif peaks == 2
    for set = 1:width(setNames)
        uselist = false(height(uptemp),1);
        nearlist = uselist;
        farlist = uselist;
        for i = 1:height(uptemp)
            if strcmp(uptemp{i,2},setNames{set})
                if uptemp{i,5} > CODlim  % cull based off of COD
                    if uptemp{i,13} > mixlim && uptemp{i,21} > mixlim
                        uselist(i,1) = 1;
                    end
                end
            end
        end
        nstat{1,set} = cell2mat(uptemp(uselist,stat1(1)));
        fstat{1,set} = cell2mat(uptemp(uselist,stat1(2)));
        output{1,set} = nstat{1,set};
        output{2,set} = fstat{1,set};
    end
    if visualize == 1
        figure()
        subplot(1,2,1)
            [nearmu,~,~,nearsem] = LAMA_CSR_err(nstat);
            bar(1:width(setNames),nearmu) , hold on 
            er = errorbar(1:width(setNames),nearmu,nearsem,nearsem);
            er.Color = [0 0 0];                            
            er.LineStyle = 'none'; 
            hold off
            ylim([0 1000])
            xticklabels(extractBefore(setNames,'_NND'))
            ylabel(strcat("Near ",stat," (nm)"))
        subplot(1,2,2)
            [farmu,~,~,farsem] = LAMA_CSR_err(fstat);
            bar(1:width(setNames),farmu) , hold on 
            er = errorbar(1:width(setNames),farmu,farsem,farsem);
            er.Color = [0 0 0];                            
            er.LineStyle = 'none'; 
            hold off
            xticklabels(extractBefore(setNames,'_NND'))
            ylim([0 1000])
            ylabel(strcat("Far ",stat," (nm)"))
    else
        [nearmu,~,~,nearsem] = LAMA_CSR_err(nstat);
        [farmu,~,~,farsem] = LAMA_CSR_err(fstat);
    end
    mu = [nearmu;farmu];
    sem = [nearsem;farsem];
end
%
end