function [img,mask] = UIA_imloader_ssv2(filename,mk,i,filenum)
    %% Function used in UIA which reads a file image, normalizes the intensity and allows for a click generated mask 
    ig = imread(filename);              % read the image file
    if size(size(ig),2)>2               % check the dimensions of the image 
        i0 = double(rgb2gray(ig));      % if there are more than 2 flatten 
    else
        i0 = double(ig);                % if not make sure pixel values are doubles anyways
    end
    qnot = i0; 
    qnot =  qnot(:); 
    qnot = (qnot-min(qnot))/(max(qnot)-min(qnot));  % normalize pixel intensity to 0,1;
    img = reshape(qnot,size(i0,1),size(i0,2));      % and reshape the image
if mk==1                                            % if we want to mask the image
        figure('Name',strcat("Image ",num2str(i)," of ",num2str(filenum)))  
        imagesc(img), hold on                       % represent the image
        title("Double Click with Mouse to Draw Area, then Press Any Key when Finished")
        tdot = 0;                                   % initialize a button press
        p = [];                                     % variable to hold mask polygon vertices
        while tdot == 0                             % as long as the user hasnt pressed a button
            [x,y] = ginput(1);                      % requerst user click input
            p = [p;[x,y]];                          
            plot(x,y,'m.','MarkerSize',14)          
            tdot = waitforbuttonpress;              
        end         
        p = round(p);                               
        p = [p;p(1,:)];                             
        [MX,MY] = meshgrid(1:size(img,2),1:size(img,1));     % generate a logical representation of the image
        [mask,~] = inpolygon(MX,MY,p(:,1),p(:,2));           % define the mask as all mesh positions inside of the polygon
        close all
    else
        mask = [];
    end
end


