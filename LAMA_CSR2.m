%% LAMA_CSR2  
%
% Created by Dylan Greening for use in RWTurner Lab UCALGARY
% 
% A function within LAMA_CSR 
% 
% which performs the randomization of
% centroid positions within the mask
%
% Input Arg:
% binaryImages =  a binary mask which is in px format. used to assess the
%   validity of randomly choosen centroid positions 
% bi_im_limits = minX,maxX,minY,maxY mask positions in px. used to speed up the randomization process 
% num_of_runs = how many randomization datasets to generate
% nm_to_px = a ratio which acts as a easy converter between nanometer and
%   pixel distances
%
% Output Arg:
% csr_cent =  a array where each cell row references a different
%   image/centroids data set and each column references a different
%   randomization run. A cell contains a n-by-2 matrix which is coordinates
%   (x,y) centroids which have been randomly placed
%%
function csr_cent = LAMA_CSR2(binaryImages,bi_im_limits,sum_of_centroids,num_of_runs,nm_to_px)
%
%sum_of_centroids = sum_of_centroids647;
%
   csr_cent = {}; %% each row is a different image, each column is a different run of randomization
    for img = 1:height(bi_im_limits)
        for run = 1:num_of_runs
            disp(strcat("Processing Image # ",num2str(img)," of ",num2str(height(bi_im_limits))," run # ",num2str(run)," of ",num2str(num_of_runs)))
            randomized_centroids = []; % generate an empty variable which we will fill for each run
            while height(randomized_centroids) < sum_of_centroids(img)
                diff_to_fill = sum_of_centroids(img) - height(randomized_centroids);
                randX = bi_im_limits(img,1) + (bi_im_limits(img,2)-bi_im_limits(img,1)) .* rand(diff_to_fill,1); % use the limits to reduce run time
                randY = bi_im_limits(img,3) + (bi_im_limits(img,4)-bi_im_limits(img,3)) .* rand(diff_to_fill,1); 
                centroids_to_check = [randX,randY]; % make the random positions look like (x,y) arrangment again
                rounded_centroids_to_check = round(centroids_to_check); % round the values so that we can check them against the mask
                approved = zeros(diff_to_fill,1);
                for i = 1:diff_to_fill
                    if binaryImages{img}(rounded_centroids_to_check(i,2),rounded_centroids_to_check(i,1)) == 1
                        approved(i,1) = 1;    
                    end
                end
                if sum(approved) > 0
                    randomized_centroids = [randomized_centroids;centroids_to_check(logical(approved),:)]; 
                end
            end
            csr_cent{img,run} = randomized_centroids./nm_to_px;
        end
    end
    %%
    %img = 1;
    %figure()
    %subplot(1,2,1)
    %imshow(binaryImages{img}), hold on
    %scatter(mca_positions{img,1}(:,3),mca_positions{img,1}(:,4),'.g')
    %scatter(mca_positions{img,2}(:,3),mca_positions{img,2}(:,4),'.b')
    %xline(bi_im_limits(img,1),'r')
    %xline(bi_im_limits(img,2),'r')
    %yline(bi_im_limits(img,3),'r')
    %yline(bi_im_limits(img,4),'r'), hold off
    %subplot(1,2,2)
    %imshow(binaryImages{img,1}), hold on
    %scatter(mca_culled647{img}(:,3) ,mca_culled647{img}(:,4) ,'.g')
    %scatter(mca_culled568{img}(:,3) ,mca_culled568{img}(:,4) ,'.b')
    %scatter(csr_cent{img,1}(:,1).*nm_to_px,csr_cent{img,1}(:,2).*nm_to_px,'.r')
    %xline(bi_im_limits(img,1),'r')
    %xline(bi_im_limits(img,2),'r')
    %yline(bi_im_limits(img,3),'r')
    %yline(bi_im_limits(img,4),'r'), hold off
    %%
end