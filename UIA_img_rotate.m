%% UIA_img_rotate
% a function used in runRefineOnOldWorkspaces and nested within UIA_NNMF.m which takes image and centroid information and 
% sythesizes 600x600 px centroid images where the Nearest neighbor is
% rotated to the left (180 degrees)
% Concieved by Wilten Nicola and adapted by Dylan Greening
function [store,rot_store,rot_store_raw,rotcent] = UIA_img_rotate(img,centroid,thresh,px,py)
    nx = size(img,1);
    ny = size(img,2);
    imgk = double(img>thresh);
    centroid = round(centroid);
    centroid(centroid(:,1)<px,:)=[]; % if the centroid data falls out of range of interest cull it
    centroid(centroid(:,2)<py,:)=[];
    centroid(centroid(:,1)+px>nx,:)=[];
    centroid(centroid(:,2)+py>ny,:)=[];
    %%% create the initial storage matrix of patches.
    nc = size(centroid,1);
    store = zeros(nc,2*px*2*py); 
    rot_store =store;  
    rand_store=store;
    storeraw = store;
    rot_store_raw = store;
    for j = 1:nc 
        cx = centroid(j,2);
        cy = centroid(j,1);
        patch = imgk(cx-px+1:cx+px,cy-py+1:cy+py);
        patchraw = img(cx-px+1:cx+px,cy-py+1:cy+py);
        store(j,:) = reshape(patch,[1,2*px*2*py]);
        storeraw(j,:) = reshape(patchraw,[1,2*px*2*py]);
    end
    md = squareform(pdist(centroid)); % distance to neighbour calcs
    md = md + 10^10*eye(nc);
    for j = 1:nc 
        [ax,ux] = sort(md(j,:));
        %% rotate the images 
        ix = ux(1);
        slope(j,1) = (centroid(j,2)-centroid(ix,2))/(centroid(j,1)-centroid(ix,1));
        angle(j,1) = -180*(atan2(centroid(j,2)-centroid(ix,2),(centroid(j,1)-centroid(ix,1))))/pi;
        rot_store(j,:) = reshape(imrotate(reshape(store(j,:),[2*px,2*py]),-angle(j,1),'crop'),[1,2*px*2*py]);
        rot_store_raw(j,:) = reshape(imrotate(reshape(storeraw(j,:),[2*px,2*py]),-angle(j,1),'crop'),[1,2*px*2*py]);
        rand_store(j,:) = reshape(imrotate(reshape(store(j,:),[2*px,2*py]),360*rand,'crop'),[1,2*px*2*py]);
        %% rotate the centroids
        i0 = find(ax<min(px,py)); 
        cent_j = [];
        cent_j = [centroid(j,:);centroid(ux(i0),:)]; 
        cent_j = [cent_j(:,1)-cent_j(1,1),cent_j(:,2)-cent_j(1,2)];
        phi = -angle(j,1)*pi/180;
        rot_j = [cos(phi),-sin(phi);sin(phi),cos(phi)];
        cent_j = cent_j*rot_j;
        rotcent{j,1} = cent_j;
    end
end