%% gammamixturemaker
%
% Created by Dylan Greening for use in RWTurner Lab UCALGARY
%
close all
clear 
clc
%
parent = "../NND-Giri-Gamma";
diroutput = dir(strcat(parent,"/"));                   % look at all of the folders in the current directory
diroutput = struct2table(diroutput);            
pullout = diroutput.isdir == 1 ;        % only look at directories
dirf = diroutput.name(pullout);         %
dirf = dirf(~ismember(dirf,'.'));       % remove unneeded names from list
dirf = dirf(~ismember(dirf,'..'));
filecount = 1; % this is a counter for the table index
% generate a master table to hold this info
% dataset / channel-pair / images / centroids / COD / mode / threshold /
% shape / rate / height / mixing portion
mastertable = table('Size',[height(dirf)*10,22],'VariableTypes', ...
    {'cellstr','cellstr','double','double','double', ...
    'double','double','double','double','double', ...
    'double','double','double','double','double', ...
    'double','double','double','double','double','double','double'},'VariableNames',...
    {'folder','channel pair','numImages','numCentroids','COD', ...
    'modeA','medianA','meanA','thresholdA','shapeA','rateA','heightA','mixingportionA', ...
    'modeB','medianB','meanB','thresholdB','shapeB','rateB','heightB','mixingportionB','Exit Flag'});
lb = [0,1,1,10,200,1,1,10]; %  phase, alpha, beta, high, phase, alpha, beta, high
ub = [300,30,400,1e7,800,15,350,1e6]; %  phase, alpha, beta, high, phase, alpha, beta, high
options = optimoptions('lsqnonlin','MaxIter',4000,'TolFun', 1e-12, 'MaxFunEvals', 2000 ,'TolX', 1e-12); 
% use a for loop to scan each list
for folder = 1:height(dirf)
    files = dir(strcat(parent,'/',dirf{folder}));
    files = {files.name};
    pullout = false(1,width(files));
    for w = 1:width(files)
        if contains(files{w},'.csv')
            pullout(1,w) = 1; 
        end
    end
    files = files(1,pullout);
    disp(sum(pullout))
    dataset = dirf{folder};
    for w = 1:width(files)
        data = readmatrix(strcat(parent,'/',dataset,'/',files{1,w}));
        images = width(data);
        pooleddata = reshape(data,[],1);
        pooleddata = pooleddata(~isnan(pooleddata));
        pooleddata = sort(pooleddata);
        centroids = height(pooleddata);                     
        h1 = histogram(pooleddata,'BinEdges',0:25:2000);
        x = h1.BinEdges;
        x = x(1:(end-1))+((x(2)-x(1))/2);
        yorig = h1.Values;
        ch_pair = extractBefore(files{w},'.csv');
        if length(ch_pair) == 2
            if strcmp(ch_pair(1),ch_pair(2))
                peak2 = 1;
            else
                peak2 = 0;
            end
        else
            if strcmp(ch_pair(1),ch_pair(2)) && ~strcmp(ch_pair(3:5),'csr')
                peak2 = 1;
            else
                peak2 = 0;
            end
        end
        gaMM1 = @(z) pdf("Gamma",x-z(1),z(2),z(3))*z(4);
        if peak2 == 1
        % set initial testing parameters
        beta = [10,20];
        alpha = [15,20];
        phase = [5,150];
        high = [1000*283,1000*283];
        gaMM2 = @(z) (pdf("Gamma",x-z(1),z(2),z(3))*z(4) + pdf("Gamma",x-z(5),z(6),z(7))*z(8));
        initial = [phase(1),alpha(1),beta(1),high(1),phase(2),alpha(2),beta(2),high(2)];
        %  phase, alpha, beta, high
       
        [solcoeff,normres,res,e_flag,~,lambda] = lsqnonlin(@(z) (gaMM2(z)-yorig),initial,lb,ub,options);
        yout = gaMM2(solcoeff);
        R2 = gfit2(yout,yorig,'8');
        gamAarea = trapz(x,gaMM1(solcoeff(1:4)));
        gamBarea = trapz(x,gaMM1(solcoeff(5:8)));
        totarea = trapz(x,yout);
        Aout = spline(x,gaMM1(solcoeff(1:4)),1:2000);
        [~,AI] = max(Aout);
        Atable = tableconverter(1:2000,Aout);
        meanA = mean(Atable);
        medianA = median(Atable);
        Bout = spline(x,gaMM1(solcoeff(5:8)),1:2000);
        [~,BI] = max(Bout);
        Btable = tableconverter(1:2000,Bout);
        meanB = mean(Btable);
        medianB = median(Btable);
        figure('Name',strcat(dataset," ",extractBefore(files{1,w},'.csv')))
        subplot(1,3,1)
        yerr = yorig-yout;
        histogram(yerr)
        [h,p] = kstest(yerr);
        title(strcat("error plot with a kstest p-value of ",num2str(round(p,4))))
        subplot(1,3,2)
        scatter(x,res), hold on
        yline(0,'-r'), hold off
        title('error along fit')
        subplot(1,3,3)
        histogram(pooleddata,'BinEdges',0:25:2000), hold on
        text(1000,max(yout),strcat("COD = ", num2str(round(R2,4))))
        text(1000,max(yout)*0.9,strcat("gammaA has a mixing portion of ",num2str(round(gamAarea/totarea,4))))
        text(1000,max(yout)*0.8,strcat("mode = ",num2str(AI)," medianA = ",num2str(round(medianA))," meanA = ",num2str(round(meanA))))
        text(1000,max(yout)*0.7,strcat("gammaB has a mixing portion of ",num2str(round(gamBarea/totarea,4))))
        text(1000,max(yout)*0.6,strcat("mode = ",num2str(BI)," medianB = ",num2str(round(medianB))," meanB = ",num2str(round(meanB))))
        plot(x,gaMM1(solcoeff(1:4)),'r-')
        plot(x,gaMM1(solcoeff(5:8)),'g-')
        plot(x,yout,'k','LineWidth',2), hold off
        set(gcf,'Units','Inch','Position',[4 4 20 5])
        saveas(gcf,strcat(parent,'/',dataset,'/',extractBefore(files{1,w},'.csv'),'.pdf'))
        %
        else
        beta = [10,20];
        alpha = [15,20];
        phase = [5,150];
        high = [1000*283,1000*283];
        initial = [phase(1),alpha(1),beta(1),high(1)];
        [solcoeff,normres,res,e_flag,~,lambda] = lsqnonlin(@(z) (gaMM1(z)-yorig),initial,lb(1,1:4),ub(1,1:4),options);
        yout = gaMM1(solcoeff);
        R2 = gfit2(yout,yorig,'8');
        figure('Name',strcat(dataset," ",extractBefore(files{1,w},'.csv')))
        subplot(1,3,1)
        yerr = yorig-yout;
        histogram(yerr)
        [h,p] = kstest(yerr);
        title(strcat("error plot with a kstest p-value of ",num2str(round(p,4))))
        subplot(1,3,2)
        scatter(x,res), hold on
        yline(0,'-r'), hold off
        title('error along fit')
        subplot(1,3,3)
        Xout = spline(x,yout,0:1:2000);
        [~,XI] = max(Xout);
        Xtable = tableconverter(x,gaMM1(solcoeff(1:4)));
        meanX = mean(Xtable);
        medianX = median(Xtable);
        fiterror = ((medianX-median(pooleddata))/median(pooleddata))*100;
        histogram(pooleddata,'BinEdges',0:25:2000),hold on
        text(1000,max(yout)*0.8,strcat("COD = ", num2str(round(R2,4))))
        text(1000,max(yout)*0.6,strcat(" mode = ",num2str(XI)," median = ",num2str(round(medianX))," mean = ",num2str(round(meanX))))
        plot(x,yout,'k','LineWidth',2), hold off
    
        set(gcf,'Units','Inch','Position',[4 4 20 5])
        saveas(gcf,strcat(parent,'/',dataset,'/',extractBefore(files{1,w},'.csv'),'.pdf'))
        end
         % generate the table with the info
        % dataset / channel-pair / images / centroids / COD / mode / threshold /
        % shape / rate / height / mixing portion
        mastertable{filecount,1} = {dataset};
        mastertable{filecount,2} = {ch_pair};
        mastertable{filecount,3} = images;
        mastertable{filecount,4} = centroids;
        mastertable{filecount,5} = R2;
        if peak2 == 0
        mastertable{filecount,6:13} = [XI,round(medianX),round(meanX),-1*solcoeff(1),solcoeff(2:4),1];
        mastertable{filecount,14:22} = [NaN,NaN,NaN,NaN,NaN,NaN,NaN,fiterror,e_flag];
        elseif peak2 ==1
        mastertable{filecount,6:22} = [AI,round(medianA),round(meanA),-1*solcoeff(1),solcoeff(2:4),gamAarea/totarea,...
            BI,round(medianB),round(meanB),-1*solcoeff(5),solcoeff(6:8),gamBarea/totarea,e_flag];
        end
        filecount = filecount + 1;
    end
end
%
writetable(mastertable,strcat(parent,"/masterNNDtableGamma.csv"));
%%
sum(mastertable{:,5})