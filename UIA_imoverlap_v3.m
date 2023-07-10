% UIA_imoverlap_v3
% function for use in UIA3
% created by Wilten Nicola 
% adapted by Dylan Greening for use in RW Turner Lab UCALGARY
function [overlap] = UIA_imoverlap_v3(files)
    nimage = size(files,2);
    for i = 1:nimage % for each imaged loaded into overlap
        img(:,:,i) = files{i}; % read the image files in a grayscale as part of a matrix
        img(:,:,i) = (img(:,:,i) - min(min(img(:,:,i))))/(max(max(img(:,:,i)))-min(min(img(:,:,i)))); %normalize the data between 1,0
    end
    overlap = ones(size(img,1),size(img,2)); % make an empty matrix of ones the same size as the images
    for j = 1:nimage
        overlap = squeeze(img(:,:,j)).*overlap; % if there is value to a index position that remains valid for the other images
    end
    overlap = nthroot(overlap,nimage); % take the nroot of the image products 
end

