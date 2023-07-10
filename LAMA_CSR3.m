%% LAMA_CSR3 
% 
% Created by Dylan Greening for use in RWTurner Lab at UCALGARY
% 
% A function within LAMA_CSR.m
%
% Takes centroid positions of real and randomized data and conducts 
% NND
% ripleysK
% radial distribution function plot
%
% Requiers functions: 
%   LAMA_CSR_pdist1
%   LAMA_CSR_nnd1
%   LAMA_CSR_nan
%   LAMA_CSR_distrange
%   LAMA_CSR_ann
%   LAMA_CSR_err
%
% Input Arg:
% mca_culled647 = an array were each cell is a n-by-2 matrix referencing
%   (x,y) coordinates of centroids within the binary mask of the soma
% csr_647 = like mca_culled but random_positions within the mask
% distrange = edges used for ripleys and radial distribution calculations
% visualize = logical if visualizations are desired
% lil_nam = two character name identifying the dataset and the labels of
%   the channels (first channel = first letter = 647)
%
% Output Arg:
%
% nnd_singles = array where each column is the nnd matrix for 
% h = figures generated to be saved in the parent script
%
%%
function [nnd_singles,h] = LAMA_CSR3(mca_culled647,csr647,mca_culled568,csr568,distrange,visualize,lil_nam)
%% evaluate point distances
% csr647 = csr_cent_647;
% csr568 = csr_cent_568;
% distrange = 0:10:2000;
% initialize some variables 
cent647 = {}; % one for centroid positions of the 647 channel
cent568 = {}; % one for centroid positions of the 568 channel
    for i = 1:height(mca_culled647) 
    cent647{i,1} = mca_culled647{i,1}(:,1:2); 
    end
    for i = 1:height(mca_culled568)
    cent568{i,1} = mca_culled568{i,1}(:,1:2);
    end   
    num_of_runs = width(csr647);
    pd_cent_647 = LAMA_CSR_pdist1(cent647); % for real centroids
    pd_cent_568 = LAMA_CSR_pdist1(cent568);
    pd_csr_647 = LAMA_CSR_pdist1(csr647);   % and randomized
    pd_csr_568 = LAMA_CSR_pdist1(csr568);
 % cross channel pds
    pd_cent_X = {}; % euclidian distances of data points
    for i = 1:height(cent647)
        for j = 1:height(cent568)
            pd_cent_X{i,j} = pdist2(cent647{i},cent568{j});
        end
    end
    pd_csr_X = {}; % euclidian distances of spatialy randomized points
    for i = 1:height(csr647)
        for j = 1:width(csr647)
            pd_csr_X{i,j} = pdist2(csr647{i,j},csr568{i,j});
        end
    end
    % this will compare the original data against every run of the random data
    pd_647orig_568csr = {};
    for i = 1:height(cent647)                               
        for j = 1:width(csr568)
            pd_647orig_568csr{i,j} = pdist2(cent647{i},csr568{i,j});
        end
    end
    pd_568orig_647csr = {};
    for i = 1:height(cent568)
        for j = 1:width(csr647)
            pd_568orig_647csr{i,j} = pdist2(cent568{i},csr647{i,j});
        end
    end
% evaluate NND with a bias to more neighbours
nnd_cent_647 = LAMA_CSR_nnd1(pd_cent_647,0); % 1D is simple and will use a like function
nnd_cent_568 = LAMA_CSR_nnd1(pd_cent_568,0);
nnd_csr_647 = LAMA_CSR_nnd1(pd_csr_647,0);
nnd_csr_568 = LAMA_CSR_nnd1(pd_csr_568,0);
%
% 2D NND
% pd_cent_X **** note only the diagonal of the array is from the same images
nnd_cent_X = {}; 
for i = 1:height(pd_cent_X) 
    for j = 1:width(pd_cent_X)
        sizeans = size(pd_cent_X{i,j});
        [~,I] = max(sizeans);
        if I == 1
            temp = sort(pd_cent_X{i,j},2);
            nnd_cent_X{i,j} = temp(:,1);
        elseif I == 2
            temp = sort(pd_cent_X{i,j},1);
            nnd_cent_X{i,j} = temp(1,:)';
        end
    end
end
% reorganize the data so only the matched images are used 
w = 1:width(nnd_cent_X);
nnd_cent_X_match = {};
for i = w
    nnd_cent_X_match{i,1} = nnd_cent_X{i,i}; % only data sets on the diagonal
end
% % data that is only not on the diagonal
nnd_cent_X_mismatch = [];
for i = 1: height(nnd_cent_X)
    for j = 1:width(nnd_cent_X)
        if i ~= j 
            nnd_cent_X_mismatch = [nnd_cent_X_mismatch;nnd_cent_X{i,j}];
        end
    end
end
% Do the CSR cross channel now 
nnd_csr_X = {}; 
for i = 1:height(pd_csr_X) 
    for j = 1:width(pd_csr_X)
        sizeans = size(pd_csr_X{i,j});
        [~,I] = max(sizeans);
        if I == 1
            temp = sort(pd_csr_X{i,j},2);
            nnd_csr_X{i,j} = temp(:,1);
        elseif I == 2
            temp = sort(pd_csr_X{i,j},1);
            nnd_csr_X{i,j} = temp(1,:)';
        end
    end
end
% now the orig-rand comparisons
nnd_568orig_647csr = {};
for i = 1:height(pd_568orig_647csr)
    for j = 1:num_of_runs
        sizeans = size(pd_568orig_647csr{i,j});
        [~,I] = max(sizeans);
        if I == 1
            temp = sort(pd_568orig_647csr{i,j},2);
            nnd_568orig_647csr{i,j} = temp(:,1);
        elseif I == 2
            temp = sort(pd_568orig_647csr{i,j},1);
            nnd_568orig_647csr{i,j} = temp(1,:)';
        end
    end
end
%
nnd_647orig_568csr = {};
for i = 1:height(pd_647orig_568csr)
    for j = 1:num_of_runs
        sizeans = size(pd_647orig_568csr{i,j});
        [~,I] = max(sizeans);
        if I == 1
            temp = sort(pd_647orig_568csr{i,j},2);
            nnd_647orig_568csr{i,j} = temp(:,1);
        elseif I == 2
            temp = sort(pd_647orig_568csr{i,j},1);
            nnd_647orig_568csr{i,j} = temp(1,:)';
        end
    end
end
% this sets up indicies for the 
[col,row] = ndgrid(1:num_of_runs,1:height(nnd_cent_647));
superindex = [reshape(row,[],1),reshape(col,[],1)];
% reformat the nnd into tables for interpretation by later scripts
nan_647_cent = LAMA_CSR_nan(nnd_cent_647,[(1:height(nnd_cent_647))',(repelem(1,height(nnd_cent_647)))']); % look at all rows of the only column
nan_568_cent = LAMA_CSR_nan(nnd_cent_568,[(1:height(nnd_cent_568))',(repelem(1,height(nnd_cent_568)))']); % look at all rows of the only column
nan_647_csr = LAMA_CSR_nan(nnd_csr_647,superindex); % look at all rows & columns
nan_568_csr = LAMA_CSR_nan(nnd_csr_568,superindex); % look at all rows & columns
nan_cent_X = LAMA_CSR_nan(nnd_cent_X_match,[(1:height(nnd_cent_X_match))',(repelem(1,height(nnd_cent_X_match)))']);
nan_csr_X = LAMA_CSR_nan(nnd_csr_X,superindex);
nnd_singles = {nan_647_cent,nan_568_cent,nan_647_csr,nan_568_csr,nan_cent_X,nan_csr_X };
%% Visualization
if visualize == 1
%% caluclate radial distributions
disp("previsualization")
%% evaluate ANN
ann_cent647 = LAMA_CSR_ann(nnd_cent_647);
ann_cent568 = LAMA_CSR_ann(nnd_cent_568);
ann_csr647 = LAMA_CSR_ann(nnd_csr_647);
ann_csr568 = LAMA_CSR_ann(nnd_csr_568);
ann_centX = LAMA_CSR_ann(nnd_cent_X);
ann_centX_trad = [] ;
ann_centX_miss = [] ;
for i = 1:height(ann_centX)
    for j = 1:width(ann_centX)
        if i == j 
        ann_centX_trad = [ann_centX_trad,ann_centX(i,j)];
        else
        ann_centX_miss = [ann_centX_trad,ann_centX(i,j)];
        end
    end
end
ann_647orig_568csr = LAMA_CSR_ann(nnd_647orig_568csr);
ann_568orig_647csr = LAMA_CSR_ann(nnd_568orig_647csr);
pooled647_568 = [];
pooled568_647 = [];
for i = 1:height(nnd_568orig_647csr)
    for j =1:width(nnd_568orig_647csr)
        pooled647_568 = [pooled647_568;nnd_647orig_568csr{i,j}];
        pooled568_647 = [pooled568_647;nnd_568orig_647csr{i,j}];
    end
end
%%
pd_cent_match = {};
for i = 1:height(pd_cent_X)
    for j = 1:width(pd_cent_X)
        if i==j
            pd_cent_match{i,1} = pd_cent_X{i,j};
        end
    end
end




    rd647pool = zeros(height(pd_cent_647),width(distrange)-1);
rd568pool = zeros(height(pd_cent_647),width(distrange)-1);
rd647ripley = rd647pool;
rd568ripley = rd568pool;
rdX = zeros(height(pd_cent_647),width(distrange)-1);
%rdXmismatch = zeros(height(pd_cent_647),width(distrange)-1);
%rdX647568 = zeros(height(pd_cent_647),width(distrange)-1);
%rdX568647 = zeros(height(pd_cent_647),width(distrange)-1);
%rdcsr647pool = zeros(num_rand*height(pd_csr_647),width(distrange)-1);
%rdcsr568pool = zeros(num_rand*height(pd_csr_647),width(distrange)-1);
%rdcsr568ripley = rdcsr647pool;
%rdcsr647ripley = rdcsr568pool;
%
for i = 1:height(pd_cent_647)
    %rd647{i,1}= LAMA_CSR_distrange(pd_cent_647{i,1},1,distrange);
    %rd568{i,1}= LAMA_CSR_distrange(pd_cent_568{i,1},1,distrange);
    rdX{i,1} = LAMA_CSR_distrange(pd_cent_match,1,distrange);
    rd647pool(i,:) = rd647{i,1}(1,:);
    rd568pool(i,:) = rd568{i,1}(1,:);
    rd647ripley(i,:) = rd647{i,1}(3,:);
    rd568ripley(i,:) = rd568{i,1}(3,:);
end
%for i = 1:height(pd_cent_X)
%    for j = 1:width(pd_cent_X)
%        rdX{j,i}= LAMA_CSR_distrange(pd_cent_X{j,i},2,distrange);
%    end
%end
%count = 1;
%for i = 1:num_rand
%   for j = 1:height(pd_csr_647)
%        rd647csr{j,i}= LAMA_CSR_distrange(pd_csr_647{j,i},2,distrange);
%        rd568csr{j,i}= LAMA_CSR_distrange(pd_csr_568{j,i},2,distrange);
%        rd647orig568csr{j,i} = LAMA_CSR_distrange(pd_647orig_568csr{j,i},2,distrange);
%        rd568orig647csr{j,i} = LAMA_CSR_distrange(pd_568orig_647csr{j,i},2,distrange);
%        rdcsr647pool(count,:) = rd647csr{j,i}(1,:);
%        rdcsr568pool(count,:) = rd568csr{j,i}(1,:);
%        rdcsr647ripley(count,:) = rd647csr{j,i}(3,:);
%        rdcsr568ripley(count,:) = rd568csr{j,i}(3,:);
%        count = count + 1;
%   end
%end
    
    
    
    
    h(1) = figure('Name',"Average Nearest Neighbour");
    temp = {ann_cent647,ann_cent568,ann_csr647(:,1),ann_csr568(:,1),ann_centX_trad,ann_centX_miss};
    [mu,~,~,sem] = LAMA_CSR_err(temp);
    bar(1:6,mu), hold on 
    er = errorbar(1:6,mu,sem.*-1,sem);
    er.Color = [0 0 0];                            
    er.LineStyle = 'none'; 
    hold off
    xticklabels({'647','568','647-CSR','568-CSR',"X-Ch","missaligned"})
    ylabel("Mean ANN (nm)")
    % 
     h(2) =figure('Name',"NND Histograms");
    subplot(2,3,1)
    histogram(nnd_cent647pool,'BinEdges',0:25:2000)
    title(strcat("Original 647 (",lil_nam(1),")"))
    subplot(2,3,2)
    histogram(nnd_cent568pool,'BinEdges',0:25:2000)
    title(strcat("Original 568 (",lil_nam(2),")"))
    subplot(2,3,3)
    histogram(nnd_cent_X_match,'BinEdges',0:25:2000)
    text(80,1600,strcat("median is ",num2str(median(nnd_cent_X_match))))
    title("cross channel")
    subplot(2,3,4)
    histogram(nnd_csr647pool,'BinEdges',0:25:2000)
    title("CSR 647")
    subplot(2,3,5)
    histogram(nnd_csr568pool,'BinEdges',0:25:2000)
    title("CSR 568")
    subplot(2,3,6)
    histogram(nnd_mismatch,'BinEdges',0:25:2000)
    title("mismatch")
    %
     h(3) =figure('Name', "Original vs CSR NND histograms");
    subplot(3,1,1)
    histogram(nnd_cent_X_match,'BinEdges',0:25:2000)
    title("Original 647 on Original 568")
    subplot(3,1,2)
    histogram(pooled647_568,'BinEdges',0:25:2000)
    title("Original 647 on CSR 568")
    subplot(3,1,3)
    histogram(pooled568_647,'BinEdges',0:25:2000)
    title("Original 568 on CSR 647")
    %
    % pool radial distribution data into means with confidence intervals
     h(4) = figure('Name',"Radial distribution functions");
    for i = 1:width(rd647pool)
        temp647{i} = rd647pool(:,i); 
        temp568{i} = rd568pool(:,i);
        temp647csr{i} = rdcsr647pool(:,i); 
        temp568csr{i} = rdcsr568pool(:,i);
    end
    centers = distrange(1:end-1)+((distrange(2)-distrange(1))/2);
    sn = 7;
    subplot(2,2,1)
    [mu,errhigh,errlow,~] = LAMA_CSR_err(temp647);
    disp(max(mu))
    plot(centers,mu,'r-'), hold on 
    plot(centers,errhigh,'r--')
    plot(centers,errlow,'r--')
    [mu,errhigh,errlow,~] = LAMA_CSR_err(temp568);
    disp(max(mu))
    plot(centers,mu,'g-')
    plot(centers,errhigh,'g--')
    plot(centers,errlow,'g--')
    hold off
    title("647 vs 568 raw radial dist function")
    subplot(2,2,2)
    [mu,errhigh,errlow,~] = LAMA_CSR_err(temp647csr);
    disp(max(mu))
    plot(centers,mu,'r-'), hold on 
    plot(centers,errhigh,'r--')
    plot(centers,errlow,'r--')
    [mu,errhigh,errlow,~] = LAMA_CSR_err(temp568csr);
    disp(max(mu))
    plot(centers,mu,'g-')
    plot(centers,errhigh,'g--')
    plot(centers,errlow,'g--')
    hold off
    title("647 vs 568 CSR radial dist function")
    subplot(2,2,3)
    [mu,errhigh,errlow,~] = LAMA_CSR_err(temp647);
    plot(centers,smooth(mu,sn),'r-'), hold on 
    plot(centers,smooth(errhigh,sn),'r--')
    plot(centers,smooth(errlow,sn),'r--')
    [mu,errhigh,errlow] = LAMA_CSR_err(temp647csr);
    plot(centers,smooth(mu,sn),'k-'), hold on 
    plot(centers,smooth(errhigh,sn),'k--')
    plot(centers,smooth(errlow,sn),'k--')
    hold off
    title("647 smoothed vs csr smoothed")
    subplot(2,2,4)
    [mu,errhigh,errlow,~] = LAMA_CSR_err(temp568);
    plot(centers,smooth(mu,sn),'g-'), hold on 
    plot(centers,smooth(errhigh,sn),'g--')
    plot(centers,smooth(errlow,sn),'g--')
    [mu,errhigh,errlow] = LAMA_CSR_err(temp568csr);
    plot(centers,smooth(mu,sn),'k-'), hold on 
    plot(centers,smooth(errhigh,sn),'k--')
    plot(centers,smooth(errlow,sn),'k--')
    title("568 smoothed vs csr smoothed")
    %%
     h(5) =figure('Name',"Ripleys(H) function");
    for i = 1:width(centers)
        temp647{i} = rd647ripley(:,i); 
        temp568{i} = rd568ripley(:,i);
        temp647csr{i} = rdcsr647ripley(:,i); 
        temp568csr{i} = rdcsr568ripley(:,i);
    end
    ax(1) = subplot(2,2,1);
    [mu,errhigh,errlow,~] = LAMA_CSR_err(temp647);
    disp(max(mu))
    plot(centers,mu,'r-'), hold on 
    plot(centers,errhigh,'r--')
    plot(centers,errlow,'r--')
    [mu,errhigh,errlow,~] = LAMA_CSR_err(temp568);
    disp(max(mu))
    plot(centers,mu,'g-')
    plot(centers,errhigh,'g--')
    plot(centers,errlow,'g--')
    hold off
    title("647 vs 568 raw ripleys function")
    ax(2) = subplot(2,2,2);
    [mu,errhigh,errlow,~] = LAMA_CSR_err(temp647csr);
    disp(max(mu))
    plot(centers,mu,'r-'), hold on 
    plot(centers,errhigh,'r--')
    plot(centers,errlow,'r--')
    [mu,errhigh,errlow,~] = LAMA_CSR_err(temp568csr);
    disp(max(mu))
    plot(centers,mu,'g-')
    plot(centers,errhigh,'g--')
    plot(centers,errlow,'g--')
    hold off
    title("647 vs 568 CSR Ripleys function")
    ax(3) = subplot(2,2,3);
    [mu,errhigh,errlow,~] = LAMA_CSR_err(temp647);
    plot(centers,smooth(mu,sn),'r-'), hold on 
    plot(centers,smooth(errhigh,sn),'r--')
    plot(centers,smooth(errlow,sn),'r--')
    [mu,errhigh,errlow,~] = LAMA_CSR_err(temp647csr);
    plot(centers,smooth(mu,sn),'k-'), hold on 
    plot(centers,smooth(errhigh,sn),'k--')
    plot(centers,smooth(errlow,sn),'k--')
    hold off
    title("647 smoothed vs csr smoothed")
    ax(4) =subplot(2,2,4);
    [mu,errhigh,errlow,~] = LAMA_CSR_err(temp568);
    plot(centers,smooth(mu,sn),'g-'), hold on 
    plot(centers,smooth(errhigh,sn),'g--')
    plot(centers,smooth(errlow,sn),'g--')
    [mu,errhigh,errlow,~] = LAMA_CSR_err(temp568csr);
    plot(centers,smooth(mu,sn),'k-'), hold on 
    plot(centers,smooth(errhigh,sn),'k--')
    plot(centers,smooth(errlow,sn),'k--')
    title("568 smoothed vs csr smoothed")
    linkaxes(ax)
else
    h = [];
end
end
