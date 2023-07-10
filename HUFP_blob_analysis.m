%% HUFP_blob_analysis
%
% Created by Dylan Greening for use in RWTurner Lab UCALGARY
%
% A function for use inside of HUFP_runone.m
%
% Input Arg:
% blob_count = scalar of blob features detected
% blob_feats = vector of indicies which define which features were blob
%   classified
% sm_fm = smoothed frame = the matrix where columns indicate features and
%   rows indicate smoothed data at pixel distance
% folder_name = location to save data
% base_file_name = name of the dataset
% px_sz = pixel size in nanometers
%
% Output Arg:
% generates pdf images and .csv files of blob analysis results
%%
function HUFP_blob_analysis(blob_count,blob_feats,sm_fm,folder_name,base_file_name,px_sz)
    %% Fit the bimodal features
    binwid = 20; % this is the width of the histogram for later
    boundary = 100:250;     % the range where gaussians will be assessed for bias to the nearest neighbour position
    decent_found = 0 ;      % make a counter for the number of decent blob fits detected 
    if blob_count > 0       % if there are any blobs to do this procedure on 
        blob_table = zeros([1,2]) ;% initialize the table x,y and a counter
        count = 1 ;
        for j = 1:width(blob_feats)     % look at all blob features
            feat = blob_feats(j);
            for index = boundary
                blob_table(count,1) = index ; % store their x and y info for the relevant range
                blob_table(count,2) = sm_fm(index,feat) ;
                count = count + 1 ;  
            end
        end 
        maxH = max(blob_table(1,:));
        % Fit a gaussian to each blob feature a
        blob_table2 =zeros(width(blob_feats),5) ;  % mean, ci, sd, ci, Rs
        for feat = 1:width(blob_feats) % this attempts to fit a single gaussian on the bounded data range
            try
            [g2fit, g2R] = fit(boundary',sm_fm(boundary,blob_feats(feat)), 'gauss1');
            g2C = confint(g2fit) ; % if the gaussian is fit and its position makes sense compared to the centering centroi
            blob_table2(feat,1) = (300 - g2fit.b1) * px_sz ;    % store mean as distance from center
            blob_table2(feat,2) = (g2C(2,2) - g2fit.b1) * px_sz ; % store confidence interval of Mean
            blob_table2(feat,3) = (g2fit.c1) * px_sz;               % store standard deviation
            blob_table2(feat,4) = (g2C(2,3) - g2fit.c1) * px_sz;    % store confidence interval of SD
            blob_table2(feat,5) = g2R.rsquare ;                     % store coeffieient of determination
            catch
                figure()    % in the case where the data could not be fit generate a figure of the problematic feature
                plot(1:600,sm_fm(1:600,blob_feats(feat)))
                title("i == ",num2str(feat))
                continue
            end
        end
        % dont use gaussians when the g2R.rsquare <90
        % of the mean or sd is outlier
        % or the mean is outside of the boundary
        [~,grah1,~,meanBoundL,meanBoundU] = rmoutliers(blob_table2(:,1),'median');
        [~,grah2,~,sdBoundL,sdBoundU] = rmoutliers(blob_table2(:,3),'median');
        grah = grah1 | grah2;
        blob_table3 = blob_table2(~grah,:); % outliers removed in blobtable3
        decent_blob_fits = zeros(1) ;       % prepare a fresh table to only hold scrutinzed blob features
        count = 1 ; 
        for i = 1:length(blob_table3(:,1))  % for every entry in blob_table3
            if blob_table3(i,5) >= 0.9 && blob_table3(i,3) < 400  % test for R2 greater than 0.9 and SD less than 400 nm
                if blob_table3(i,1) > (300 - max(boundary))*px_sz && blob_table3(i,1) < (300-min(boundary))*px_sz % ensure mean in within bounds
                    decent_blob_fits(count) = i ; % if all scruitny is passed store in the scrutinized table
                    count = count + 1 ; 
                end
            end
        end                                 
        if length(decent_blob_fits) ~= 1 && decent_blob_fits(1) ~= 0 % If there are still blob features to evaulate 
            decent_blobs = zeros(length(decent_blob_fits),5);        % Make a new table which only evaluates the well fitted blobs by looking at their gaussian coefficients
            row_count = 1 ;
            for i = decent_blob_fits                                
                decent_blobs(row_count,:) = blob_table3(i,:);
                row_count = row_count + 1 ; 
            end
            % make images of all of the distributions from the isolated data
            close all
            subplot(1,4,1), hold on 
            plot(blob_table2(:,1), blob_table2(:,3),'o')
            axis([-20 1600 0 400])
            ylabel('sd of fitted gaussian (nm)')
            xlabel('mean of fitted gaussian (nm)')
            xline(meanBoundL,'--r')
            xline(meanBoundU,'--r')
            yline(sdBoundL,'--r')
            yline(sdBoundU,'--r')
            hold off
            subplot(1,4,2)
            plot(decent_blobs(:,1), decent_blobs(:,3),'o')
            axis([-20 1600 0 400])
            ylabel('sd of fitted gaussian (nm)')
            xlabel('mean of fitted gaussian (nm)')
            xline(meanBoundL,'--r')
            xline(meanBoundU,'--r')
            yline(sdBoundL,'--r')
            yline(sdBoundU,'--r')
            xline((300 - max(boundary))*px_sz ,'-r','LineWidth',2)
            xline((300-min(boundary))*px_sz, '-r' ,'LineWidth',2)
            hold off
            subplot(1,4,3)
            hold on 
            for i = decent_blob_fits
                feat = blob_feats(i);
                plot(200:-1:50,sm_fm(boundary,feat), 'Color','black')
                xlabel('px displacement')
            end  
            ylim([0,0.08])
            hold off
            subplot(1,4,4)
            histA = histogram(decent_blobs(:,1),'BinWidth', binwid, 'FaceColor', 'b') ; 
            axis([-10 1500 0 (max(histA.BinCounts)*1.1)])
            ylabel("count of gaussian mean")
            xlabel("mean of fitted gaussian (nm)")
            set(gcf, "units", 'inches', 'position', [2 2 20 4])
            fig2_name = strcat(folder_name,'isolated images/',base_file_name,"_well_fitted_blobs.pdf") ; 
            saveas(gcf, fig2_name,'pdf')
            decent_found = 1 ;
        else
            decent_blobs = 0;
        end
        %% use GMM to define the position of the fitted blobs
        if decent_found == 1 && height(decent_blobs) >= 5
            satisfied = 'n' ; 
            while satisfied == 'n'
                %% 
                %prompt = "How many gaussians would you like to check for the GMM?   ";
                g_num = 1;
                if g_num > 0                                    % If more than 0 gaussian were selected
                    close all
                    decentmusd = zeros(length(decent_blobs(:,1)),2) ;  % initialize a matrix
                    decentmusd(:,1) = decent_blobs(:,1) ; 
                    decentmusd(:,2) = decent_blobs(:,3) ; 
                    for i = 1:g_num
                        gm = fitgmdist(decentmusd,i) ;                 % fit the distribution for every number of gaussians
                        if i == 1
                            first_fit_GMM = gm ; 
                        end
                        subplot(1,g_num,i)
                        scatter(decentmusd(:,1),decentmusd(:,2),40,'.') % Scatter plot with points of size 10
                        axis([-10 1500 0 400])
                        hold on
                        gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(gm,[x0 y0]),x,y);
                        fcontour(gmPDF,[-10 1500 0 400])
                    end
                    hold off
                    drawnow
                    
                end
                satisfied = 'y' ; 
                %prompt = "Are you satisfied with the fit? y or n   " ;
                %satisfied = input(prompt,'s') ; 
            end
            %prompt = "How many gaussians do you actually want to use?   ";
            g_num = 1 %input(prompt);
            close all
            subplot(1,3,1)
            scatter(decentmusd(:,1),decentmusd(:,2),40,'.') % Scatter plot with points of size 10
            axis([-10 1500 0 400])
            %
            subplot(1,3,2)
            hold on 
            histA = histogram(decent_blobs(:,1),'BinWidth', binwid, 'FaceColor', 'b') ; 
            axis([0 1500 0 (1.1*max(histA.BinCounts))])
            ylabel("count of gaussian mean")
            xlabel("mean of fitted gaussian (nm)")
            gm = fitgmdist(decentmusd(:,1),g_num) ;
            for i = 1:g_num  
            xline(gm.mu(i,1),'--g',gm.mu(i,1))
            end 
            hold off
            subplot(1,3,3)
            %
            hold on
            histB = histogram(decent_blobs(:,1),'BinWidth',binwid,'FaceColor','b') ; 
            axis([-10 1500 0 (max(histA.BinCounts)*1.1)])
            ylabel("count of gaussian mean")
            xlabel("mean of fitted gaussian (nm)")
            [decenty,decentx] = ksdensity(decent_blobs(:,1),1:1500,'Bandwidth',35) ;       % do a density plot of the histogram data
            plot(decentx,(max(histA.BinCounts)/max(decenty))*decenty)       % plot it scaled to the histogram height
            [~,di] = findpeaks((max(histA.BinCounts)/max(decenty))*decenty,'MinPeakProminence',1,'MinPeakHeight',1.5)   ;                                % find peaks in the density plot
            for i = 1:width(di)  
            xline(di(i),'-r',di(i))
            end         
            hold off          
            %
            set(gcf, "units", 'inches', 'position', [2 2 20 4])
            fig2_name = strcat(folder_name,'isolated images/',base_file_name,"_GMM_blobs.pdf") ; 
            saveas(gcf, fig2_name,'pdf')
        else
            g_num = 0;
        end
        %
        try
            writematrix(decent_blobs(:,1), strcat(folder_name,'isolated data/',base_file_name,'_origin.csv')) ; 
        catch
            writecell(decent_blobs(:,1), strcat(folder_name,'isolated data/',base_file_name,'_origin.csv')) ; 
        end
        %
        eq = zeros((1+g_num),6);
        mu1 = mean(decent_blobs(:,1));
        sd1 = std(decent_blobs(:,1));
        n1 = height(decent_blobs(:,1));
        min1 = min(decent_blobs(:,1));
        max1 = max(decent_blobs(:,1));
        sem1 = sd1/sqrt(n1);
        eq(1,:) = [n1,mu1,sd1,sem1,min1,max1] ;
        for i = 1:g_num
             eq((1+i),1) = gm.ComponentProportion(1,i);
             eq((1+i),2) = gm.mu(i,1);
             eq((1+i),3) = gm.Sigma(1,1,i);
        end
        col_names = {'n' 'mean' 'sd' 'sem' 'min' 'max'} ;
        eq = array2table(eq);
        eq.Properties.VariableNames = col_names;
        writetable(eq, strcat(folder_name,'isolated data/',base_file_name,'_stat.csv')) ; 
    else
        error("blob_count is less than or equal to zero which is incompatible with blob_analysis")
    end
   close all
   disp(eq)
end