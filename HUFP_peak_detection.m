%% HUFP_peak_detection
%
% Created by Dylan Greening for use in RWTurner Lab UCALGARY
%
% A function used in HUFP_runone.m 
%
% Used to classify NNMF outputs into ladder or blob features, 
% also generates a smoothed version of the data sets as well as indicies where peaks were found.  
% 
% Input Arg:
% x0 = pixel positions in the centoid images
% raw = raw intensity information from each position in the features
% feat_num = how many features are in the dataset
% feat_per_cell = how many features were drawn from each cell by NNMF
%
% Output Arg: 
% sm_fm = smoothed frame = a matrix where columns refer to different
%   features referencing 'raw' after a spline smoothing
% pk_fm = frame refenerencing raw and sm_fm defining the positions of peaks
% ladder_feats = vector of indicies defining ladder classified features
% blob_feats = vector of indicies defining blob classified features
% neither_feats = vector of indicies defining unclassified features
% blob_count = how many blob feats were detected
%%
function [sm_fm,pk_fm,ladder_feats,blob_feats,neither_feats,blob_count] = HUFP_peak_detection(x0,raw,feat_num,feat_per_cell)
    p_MIN = 0.001     ; % promemnence threshold
    h_MIN = 0.05      ; % height threshold
    d_MIN = 12        ; % distance from other peak threshold in pixels
    sm_fm = zeros([600,feat_num]); % initialized dataframe for smoothed data
    pk_fm = zeros([600,feat_num]); % data frame for peak positions
    % Fill sm_fm with the new data
    for i = 1:feat_num
        [~, values] = spaps(x0,raw(:,i),0.000005); 
        sm_fm(:,i) = values;
    end
    % Plot peaks to see how things are looking
    for i = 1:feat_num
        temparray = sm_fm(:,i); 
        [pks,index] = findpeaks(temparray,'MinPeakDistance',d_MIN,'MinPeakProminence',max(temparray)*p_MIN,'MinPeakHeight',(max(temparray)*h_MIN));
        % Build a table of peaks 
        for j = 1:length(index) % look at each peaks index 
            pk_fm(index(j),i) = pks(j) ; % store peak heights at the proper positions
        end
    end
    % Get user request if plots need to be made
    prompt = 'Do you want to plot all of the detected peaks on the smoothed data ? y or n       ';
    jam = input(prompt,'s') ; 
    if (jam  == 'Y') || (jam == 'y')
        n_cell = feat_num / feat_per_cell ;
        for cell = 1:n_cell
            cell_nam = num2str(cell) ;
            figure('Name',strcat('Image #',cell_nam))
            for feat = 1:feat_per_cell
                subplot(feat_per_cell,1,feat)
                count = 1 ; 
                peak_spot = zeros(1) ; 
                feature = ((cell-1) * feat_per_cell +feat) ; 
                plot(x0,sm_fm(:,feature)), hold on
                plot(x0,raw(:,feature),'r')
                %
                temp_vec = pk_fm(:,feature) > 0 ;
                if sum(temp_vec) > 0     
                    for i = x0
                        if temp_vec(i) == 1 
                        xline(i,'black')
                        end
                    end
                end
                %
                scatter(peak_spot, repelem((max(sm_fm(:,feature))*0.5),length(peak_spot)))
            end
        end   
    end
    % Find the portion of ladder features and the portion of blob features
    ladder_count  = 0       ; 
    blob_count = 0          ;
    neither_count = 0       ;
    blob_feats = zeros(1)   ;
    ladder_feats = zeros(1) ;
    neither_feats = zeros(1);
    range_to_use = 150:400 ;                    % look only at indicies in this range
    for feature = 1:feat_num                    % look at each feature
        B = pk_fm(range_to_use,feature) > 0 ;   % if the positions have a value greater than 0 there is a peak
        temp_peak_count = sum(B) ;              % add up how mant peaks were detected
        if temp_peak_count == 2 || temp_peak_count == 3 
            blob_count = blob_count + 1;        % if there are three or fewer peaks update blob counter
            blob_feats(blob_count) = feature ;
        elseif temp_peak_count >= 4
            ladder_count = ladder_count + 1;    % if there are 4 or more update ladder counter
            ladder_feats(ladder_count) = feature ; 
        else 
            neither_count = neither_count + 1;
            neither_feats(neither_count) = feature;
        end
    end
    %
   %
    ladder_portion = ladder_count / feat_num ;
    blob_portion = blob_count / feat_num  ;   
    neither_portion = neither_count / feat_num ;
    % Check classification
    prompt = 'Do you want to plot classified peaks ? y or n       ';
    jam = input(prompt,'s');
    if (jam  == 'Y') || (jam == 'y')
        close all hidden
        prompt = 'Do you want to plot master features ? y or n       ';
        jam = input(prompt,'s');
        if (jam  == 'Y') || (jam == 'y') % this plots a master feature which is the mean of each classified group
            mean_blob = zeros(600,1);
            for feat = 1:blob_count
                mean_blob = mean_blob+sm_fm(:,blob_feats(feat));
            end
            mean_blob = mean_blob/blob_count;
            mean_lad = zeros(600,1);
            for feat = 1:ladder_count
                mean_lad = mean_lad+sm_fm(:,ladder_feats(feat));
            end
            mean_lad = mean_lad/ladder_count;
            subplot(1,2,1)
            title("blob master feature")
            plot(x0,mean_blob,'-b'), hold on 
            [~ , ploc] = findpeaks(mean_blob);
            for i = 1:height(ploc)
                xline(ploc(i,1),'-r')
            end
            hold off
            subplot(1,2,2)
            title("ladder master feature")
            plot(x0,mean_lad,'-b'), hold on
            [~ , ploc] = findpeaks(mean_lad);
            for i = 1:height(ploc)
                xline(ploc(i,1),'-r')
            end
            hold off
        end
        % blob feat
        prompt = "enter 'y' to plot blob features     " ; 
        jam1 = input(prompt,'s') ; 
        if (jam1  == 'Y') || (jam1 == 'y')
            if blob_count < 30
                skip = 1;
            elseif blob_count < 50
                skip = 2;
            else
                skip = 4;
            end
            for i = 1:skip:width(blob_feats)
                feat = blob_feats(1,i);
                figure('Name','Blob feature')
                plot(x0, raw(:,feat), 'r'), hold on 
                plot(x0, sm_fm(:,feat), 'b')
                ylim([0 (1.1*max(raw(:,feat)))])
                temp_peaks = pk_fm(:,feat)>0 ;
                for index = 1:length(temp_peaks)
                    if temp_peaks(index,1) == 1
                        xline(index,'-r') 
                    end
                end
                hold off
            end
        end
        % ladder feat plotter
        prompt = "enter 'y' to plot ladder features    " ; 
        jam2 = input(prompt,'s') ; 
        if (jam2  == 'Y') || (jam2 == 'y')
            if ladder_count < 30
                skip = 1;
            elseif ladder_count < 50
                skip = 2;
            else
                skip = 4;
            end
            for i = 1:skip:width(ladder_feats)
                feat = ladder_feats(1,i);
                figure('Name','Ladder feature')
                title(strcat('Feature # ', num2str(feat)))
                plot(x0,raw(:,feat),'r',x0,sm_fm(:,feat),'b'), hold on 
                temp_peaks = pk_fm(:,feat)>0 ;
                ylim([0 (1.1*max(raw(:,feat)))])
                for index = 1:length(temp_peaks)
                    if temp_peaks(index,1) == 1
                        xline(index,'-r') 
                    end
                end
                hold off
            end
        end
 
    end
end