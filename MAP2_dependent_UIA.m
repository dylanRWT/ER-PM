mapfiles = dir('*488.tif');  
mapfiles = {mapfiles.name}; % select the map2 image file
%%
for i = 1:7
    figure()
    subplot(2,1,1)
    imagesc(img_array_red{1,i})%,hold on
    map2image = imread(mapfiles{i}); % read the map file
    map2image = imresize(map2image, 20);  
    subplot(2,1,2)
    imagesc(map2image)
end
%%
neg_array = {};
pos_array = {};
for i = 1:7
    map2image = imread(mapfiles{i}); % read the map file
    map2image = imresize(map2image, 20);  
    map2image = double(map2image)/double(max(map2image,[],'all'));
    q = sort(reshape(map2image,size(map2image,1)*size(map2image,2),1));
    w = q(round(height(q)*0.5),1);
    % sort masks
    neg_array{i} = zeros(width(map2image),height(map2image));
    for eachcol = 1:width(map2image)
        for eachrow = 1:height(map2image)
            if map2image(eachrow,eachcol) < w
                neg_array{i}(eachrow,eachcol) = img_array_red{1,i}(eachrow,eachcol);
            else
                pos_array{i}(eachrow,eachcol) = img_array_red{1,i}(eachrow,eachcol);
            end
        end
    end
    % MAP2 positive mask
    pos_array{i} = pos_array{i}.* abs(mask_array{i}-1);
end
%% check that the masking worked properly
for i = 1:7
    figure()
    subplot(1,3,1)
    map2image = imread(mapfiles{i}); % read the map file
    map2image = imresize(map2image, 20);  
    imagesc(map2image)
    subplot(1,3,2)
    imagesc(pos_array{i})
    subplot(1,3,3)
    imagesc(neg_array{i})
end
%%
P = [];
N = [];
nof = 6;
for i = 1:7 % for each file
    [~,centroidX,~] = UIA_region_props(pos_array{i},thresh);
    [~,~,rotrawX] = UIA_img_rotate(pos_array{i},centroidX,thresh,px,py);
    %do NNMF
    [~,h] = nnmf(rotrawX,nof);
    for j = 1:nof 
        a = reshape(h(j,1:(2*px*2*py)),[2*px,2*py]);
        P(:,end+1) = a(px,:)';
    end
    [~,centroidX,~] = UIA_region_props(neg_array{i},thresh);
    [~,~,rotrawX] = UIA_img_rotate(neg_array{i},centroidX,thresh,px,py);
    %do NNMF
    [~,h] = nnmf(rotrawX,nof);
    for j = 1:nof 
        a = reshape(h(j,1:(2*px*2*py)),[2*px,2*py]);
        N(:,end+1) = a(px,:)';
    end
end
%%
writematrix(P,"Spos6X.csv")
writematrix(N,"Sneg6X.csv")
%%
save("Positive&Negative_MAP2_UIA.mat")
