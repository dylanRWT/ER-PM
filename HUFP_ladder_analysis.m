%% HUFP_ladder_analysis
% a function used in HUFP_runone.m takes information about peak position
% assesses ladder classified features for interpeak distances and fits 
% them with a damped oscillatory model using the functions  
%   HUFP_osscilation_fitter_auto.m
%   HUFP_osscilation_fitter_manual.m
function [best_period_px,decent_displacements,fixed_rungs,rung_sum,initialValueTable,best_model,displacements] = HUFP_ladder_analysis(x0,pk_fm,ladder_feats,lwcp,upcp, folder_name, base_file_name, px_sz)  
    dcount = 1 ; 
    displacements = [];
    for i = 1:width(ladder_feats)  % 
        feat = ladder_feats(1,i) ; % look at each feature
        array1 = zeros(1) ;    % initialize an array to temp hold peak positions
        count = 1 ;            % initialize a counter for the temp array
        for index =x0     % look at each index of the feature
            if pk_fm(index,feat) > 0 
                array1(count) = index ; 
                count = count + 1 ;
            end
        end
        while array1 > 1 
            for k = 1:1:(length(array1) - 1) 
                displacements(dcount) = array1(length(array1)) - array1(k) ; 
                dcount = dcount + 1 ; 
            end
            array1(end) = [] ;
        end
    end
    % Convert histogram into (x,y) position
    y0 = zeros(1,300);
    for index = 1:length(y0)
        y0(index)= sum(displacements==index);
    end
    t0 = smooth(y0);    
    t1 = smooth(t0);
    [~, i1] = findpeaks(t1,'MinPeakDistance',5) ; % find the peaks of the double smoothed data 
    % whittle this down to the first three peaks 
    i1 = i1(1:3,1);
    initial_period = 0.5 * (i1(2) - i1(1) + i1(3) - i1(2));
    % Cull the data so that there is not fitting performed on the first 10 positions or final 100
    xCull = lwcp:1:upcp  ;
    y_smooth = t0(xCull)' ; % note this is transposed so things will work
    y_raw = y0(xCull)     ;
    nbin= 300;
    histogram(displacements,nbin)
    %%
    % We need resnorm to choose the best fit situation which will likely
    % occur on the smoothed data, as well as the coefficeients and the
    % model values
    % get results for raw data with manual settings
    [rmCoeff,rmJacobian,rmResnorm,rmExitFlag,rmSolout,rmInput] = HUFP_osscilation_fitter_manual(xCull,y_raw);
    % raw data auto intials
    [raCoeff,raJacobian,raResnorm,raExitFlag,raSolout] = HUFP_osscilation_fitter_auto(xCull,y_raw,initial_period);
    % smoothed manual
    [smCoeff,smJacobian,smResnorm,smExitFlag,smSolout,smInput] = HUFP_osscilation_fitter_manual(xCull,y_smooth);
    % smoothed auto intials
    [saCoeff,saJacobian,saResnorm,saExitFlag,saSolout] = HUFP_osscilation_fitter_auto(xCull,y_smooth,initial_period);
    % The best fitted model is used
    allResnorm = [rmResnorm,raResnorm,smResnorm,saResnorm];
    [~,best_fit_row] = min(allResnorm,[],"all");
    ismemMatrix = ismember(allResnorm,allResnorm(best_fit_row));
    if sum(ismemMatrix(:,1))>0
        disp("The best model fit is from raw data and manual inputs")
        [~,temp] = min(rmResnorm);
        best_coeff = rmCoeff;
        best_jacob = rmJacobian;
        best_model = [(y_raw/max(y_raw));(y_smooth/max(y_smooth));rmSolout(temp,:)] ;%
        exitflag = rmExitFlag;
        BM = 1 ;
    elseif sum(ismemMatrix(:,2))>0
        disp("The best model fit is from raw data and auto inputs")
        [~,temp] = min(raResnorm);
        best_coeff = raCoeff;
        best_jacob = raJacobian;
        best_model = [(y_raw/max(y_raw));(y_smooth/max(y_smooth));raSolout(temp,:)] ;%
        exitflag = raExitFlag;
         BM = 2;
    elseif sum(ismemMatrix(:,3))>0
        disp("The best model fit is from smooth data and manual inputs")
        [~,temp] = min(smResnorm);
        best_coeff = smCoeff;
        best_jacob = smJacobian;
        best_model = [(y_raw/max(y_raw));(y_smooth/max(y_smooth));smSolout(temp,:)]; %
        exitflag = smExitFlag;
         BM = 3;
    elseif sum(ismemMatrix(:,4))>0
        disp("The best model fit is from smooth data and auto inputs")
        [~,temp] = min(saResnorm);
        best_coeff = saCoeff;
        best_jacob = saJacobian;
        best_model = [(y_raw/max(y_raw));(y_smooth/max(y_smooth));saSolout(temp,:)] ;%
        exitflag = saExitFlag;
         BM = 4;
    else
        error("There was an error in the best model fit* HUFP_ladder_analysis")
    end 
    initialValueTable = zeros(2,6) ; % initial value table 
    initialValueTable(1,:) = rmInput;
    initialValueTable(2,:) = smInput;
    initialValueTable(3,:) = best_coeff;
    initialValueTable = array2table(initialValueTable);
    initialValueTable.Properties.VariableNames = {'Amplitude' 'Amp Decay' 'Period' 'Phase Shift' 'Baseline' 'Baseline Decay'} ;
    initialValueTable.Properties.RowNames = {'Raw Input' 'Smooth Input' 'Best Coeff'} ;
    writetable(initialValueTable, strcat('/Users/turnerlab/Desktop/ER-PM junctions/HUFP/UIA_v2/HUFP outputs/ladder data/model coefficients/',base_file_name,'_coeff.csv'),'WriteRowNames' , true) 
    best_model(4,:) = best_model(1,:) - best_model(3,:); % raw error
    best_model(5,:) = best_model(2,:) - best_model(3,:); % smooth error
    best_model(6,:) = best_model(4,:).^2 ; % squared error raw 
    best_model(7,:) = best_model(5,:).^2 ; % squared error smooth
    ssq_raw = sum(best_model(6,:));
    ssq_smooth = sum(best_model(7,:));
    MSE_raw = mean(best_model(6,:));
    MSE_smooth = mean(best_model(7,:));
    R2_raw = gfit2(best_model(1,:),best_model(3,:),'8');
    R2_smooth = gfit2(best_model(2,:),best_model(3,:),'8');
    % applied to raw data
    best_model = array2table(best_model);
    splitvars(best_model);
    best_model.Properties.RowNames = {'y_raw' ' y_smooth' 'y_model' 'error raw' 'error_smooth' 'squared error raw' ' squared error smooth'} ;
    writetable(best_model, strcat('/Users/turnerlab/Desktop/ER-PM junctions/HUFP/UIA_v2/HUFP outputs/ladder data/model fit/',base_file_name,'_',num2str(BM),'_fit.csv'),'WriteRowNames',true) ;
    residual_info = array2table([ssq_raw,MSE_raw,R2_raw;ssq_smooth,MSE_smooth,R2_smooth;ssq_raw/ssq_smooth,MSE_raw/MSE_smooth,R2_raw/R2_smooth]);
    residual_info.Properties.VariableNames = {'SSQ' 'MSQ' 'COD'};
    residual_info.Properties.RowNames = {'raw' 'smooth' 'ratio'};
    writetable(residual_info,strcat('/Users/turnerlab/Desktop/ER-PM junctions/HUFP/UIA_v2/HUFP outputs/ladder data/model fit/',base_file_name,'_',num2str(BM),'_resinfo.csv'),'WriteRowNames',true)
    adjusted_displacements = displacements(displacements < upcp).*px_sz ; 
    w = figure('Name',"Model Outputs");
    subplot(2,3,1) % raw data with the model
        plot(xCull,best_model{1,:},'-b',xCull,best_model{3,:},'-r',"LineWidth",2)
        legend('Raw Data' ,'Model')
        title("Raw with model")
        text(60,0.90,strcat(sprintf('%.2f',best_coeff(3)*px_sz)," nm"))
        text(60,0.80,sprintf('%.4f',R2_raw))
    subplot(2,3,2) % smooth data with the model
        plot(xCull,best_model{2,:},'-b',xCull,best_model{3,:},'-r',"LineWidth",2)
        legend('Smooth Data' ,'Model')
        title("Smooth with model")
        text(60,0.90,strcat(sprintf('%.2f',best_coeff(3)*px_sz)," nm"))
        text(60,0.80,sprintf('%.4f',R2_smooth))
    subplot(2,3,3) % model over the raw histogram
        nbin = 56/2;
        hF1 = histogram(adjusted_displacements,nbin,'Normalization','probability'), hold on ;
        new_mean = mean(hF1.BinCounts) / width(hF1.Data) ; 
        axis([(lwcp*px_sz) (upcp*px_sz) 0 (max(hF1.Values)*1.2)])
        xlabel('Interpeak distances (nm)')
        ylabel('Probability')
        if BM == 3 || BM == 4
            if BM == 3
            title("Model Generated from smooth data using manual inputs")
            else
            title("Model Generated from smooth data using automated inputs")
            end
            plot(xCull*px_sz,(best_model{3,:}./(mean(best_model{3,:})/new_mean)),'r',"LineWidth",2), hold off
            subplot(2,3,6)
            scatter(xCull,best_model{5,:},'xk'), hold on 
            ylabel("error data-model")
            ylim([-1 1])
            yline(0,'-r'), hold off
        elseif BM == 1 || BM == 2
            if BM == 1
            title("Model Generated from raw data using manual inputs")
            else
            title("Model Generated from raw data using automated inputs")
            end
            plot(xCull*px_sz,(best_model{3,:}./(mean(best_model{3,:})/new_mean)),'r',"LineWidth",2), hold off
            subplot(2,3,6)
            scatter(xCull,best_model{4,:},'xk'), hold on 
            ylabel("error data-model")
            ylim([-1 1])
            yline(0,'-r','LineWidth',2), hold off
        end
    subplot(2,3,4)
        h2 = histogram(best_model{4,:},'BinEdges',linspace(-1,1,50)), hold on 
        [SW,SF] = swft(best_model{4,:}');
        text(-1,max(h2.Values)*0.8,sprintf("Shaprio-Wilks p-value = %g",SW{2,7}))
        text(-1,max(h2.Values)*0.6,sprintf("Shaprio-Francia p-value = %g",SF{2,7}))
        [~, p] = ttest(best_model{4,:}');
        text(-1,max(h2.Values)*0.4,sprintf("ttest mean of zero \np-value = %.3g",p))
        xline(0,'-r','LineWidth',2), hold off
        xlim([-1.2 1.2])
        ylim([0 max(h2.Values)*1.2])
    subplot(2,3,5)
        h3 = histogram(best_model{5,:},'BinEdges',linspace(-1,1,50)), hold on 
        [SW,SF] = swft(best_model{5,:}');
        text(-1,max(h3.Values)*0.8,sprintf("Shaprio-Wilks p-value = %g",SW{2,7}))
        text(-1,max(h3.Values)*0.6,sprintf("Shaprio-Francia p-value = %g",SF{2,7}))
        [~, p] = ttest(best_model{5,:}');
        text(-1,max(h3.Values)*0.4,sprintf("ttest mean of zero \np-value = %.3g",p))
        xline(0,'-r','LineWidth',2), hold off
        xlim([-1.2 1.2])
        ylim([0 max(h3.Values)*1.2])
    set(w,'Units','Normalized','OuterPosition',[0 0 1 1])
    fig1_name = strcat('/Users/turnerlab/Desktop/ER-PM junctions/HUFP/UIA_v2/HUFP outputs/ladder images/model fit/',base_file_name,".pdf");
    print(fig1_name,'-dpdf','-bestfit')
    close all
    best_period_px = best_coeff(3);
    decent_displacements = zeros(width(ladder_feats),1);
    fixed_pos = [];
    for i = -8:1:8
        %fixed_pos(end+1,1) = round(300+(i*best_period_px)-1);
        fixed_pos(end+1,1) = round(300+(i*best_period_px));
       % fixed_pos(end+1,1) = round(300+(i*best_period_px)+1);
    end
    fixed_rungs = zeros(height(fixed_pos),width(ladder_feats));
    for feat = 1:width(ladder_feats)
        temp = pk_fm(:,ladder_feats(feat));
        for i = 1:height(fixed_pos)
            index = fixed_pos(i,1);
            if (temp(index,1) > 0 || temp(index-1,1) > 0) || temp(index+1,1) > 0 
                fixed_rungs(i,feat) = 1;
            end
        end
        pk_positions = [];
        for i = 150:450
            if temp(i) > 0
            pk_positions(end+1,1) = i;
            end
        end
        count = 0;
        for i = height(pk_positions):-1:2
            mini_disp = (pk_positions(i) - pk_positions((i-1),1));
            if (mini_disp > (round(best_period_px) - 1)) && (mini_disp < (round(best_period_px) + 1)) 
                count = count + 1;
            end
        end
        decent_displacements(feat,1) = count;
    end
    rung_sum = zeros(17,1);
    for i = 1:17
        rung_sum(i,1) = sum(fixed_rungs(i,:));
    end
end
