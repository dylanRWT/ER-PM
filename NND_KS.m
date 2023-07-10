%% Kosmolov-Smirnov Test Script for NND distributions
% NND_KS
clear all
clc
close all
% select two or more files to compare
[file,path] = uigetfile("../*.csv",'multiselect','on'); 
runcols = 0; %input("Do you want to run individual columns as a seperate fit? 1 or 0    :");
% load each data set into an array
    data = {};
    for i = 1:width(file)
        data{i,1} = readmatrix(strcat(path,file{i}));
    end
% define an insignificant P-value from the CSR data
if runcols == 1
pCSR = [];
ksCSR = [];
pdata= [];
ksdata = [];
pX = [];
ksX = [];
count = 1;
countX = 1;
for img = 1:(width(data{2})-1)
    for j = (img+1):width(data{2})
        [h,pCSR(count,1),ksCSR(count,1)] = kstest2(data{1}(:,img),data{1}(:,j));
        [~,pdata(count,1),ksdata(count,1)] = kstest2(data{2}(:,img),data{2}(:,j));
        count = count + 1;
    end
    [~,pX(countX,1),ksX(countX,1)] = kstest2(data{1}(:,img),data{2}(:,img));
    countX = countX+1;
end

figure('Name','CSR vs Data distributions')
hold on
for img = 1:width(data{2})
subplot(2,width(data{2}),img)
histogram(data{1}(:,img), 'BinEdges',0:25:2000)
subplot(2,width(data{2}),img+width(data{2}))
histogram(data{2}(:,img), 'BinEdges',0:25:2000)
end
hold off
% pooled KS
centroidcount = sum(~isnan(data{2}),'all');
if runcols ==1
pCSRlog = log10(pCSR)*-1;
pdatalog = log10(pdata)*-1;
pXlog = log10(pX)*-1;
subplot(1,2,1)
scatter(repelem(1,height(pCSR)),pCSRlog, 'jitter','on', 'jitterAmount',0.3), hold on
scatter(repelem(2,height(pdata)),pdatalog, 'jitter','on', 'jitterAmount',0.3)
scatter(repelem(3,height(pX)),pXlog, 'jitter','on', 'jitterAmount',0.3)
yline(log10(pBig)*-1,'-r'), hold off
xticks([1,2,3])
xticklabels({extractBefore(file{1},"_NND"),extractBefore(file{2},"_NND"),"SingleX"})
xlim([0,4])
title("Log-negative p-values")
subplot(1,2,2)
scatter(repelem(1,height(pCSR)),ksCSR, 'jitter','on', 'jitterAmount',0.3), hold on
scatter(repelem(2,height(pdata)),ksdata, 'jitter','on', 'jitterAmount',0.3)
scatter(repelem(3,height(pX)),ksX, 'jitter','on', 'jitterAmount',0.3)
yline(log10(ksBig)*-1,'-r'), hold off
xticks([1,2,3])
xticklabels({extractBefore(file{1},"_NND"),extractBefore(file{2},"_NND"),"SingleX"})
xlim([0,4])
title("ks-stats")
end
else
  A =  reshape(data{1},[],1);
  A = A(A<250);
  B =  reshape(data{2},[],1);
  B = B(B<250);
figure()
subplot(2,2,1)
histogram( reshape(data{1},[],1),'BinEdges',0:25:2000)
subplot(2,2,2)
histogram( A,'BinEdges',0:25:250)
subplot(2,2,3)
histogram( reshape(data{2},[],1),'BinEdges',0:25:2000)
subplot(2,2,4)
histogram( B,'BinEdges',0:25:250)
[~,pBig,ksBig] = kstest2(A,B);
[pW,~] = ranksum(A,B);
disp(file{1})
disp(file{2})
disp(num2str(log10(pBig)*-1))
disp(num2str(ksBig))
disp(num2str(log10(pW)*-1))
end