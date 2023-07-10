%% LAMA_CSR_distrange
%
% Created by Dylan Greening for use in RWTurner Lab UCALGARY
%
% A function used in LAMA_CSR.m - > LAMA_CSR3.m
%
% takes pdist information and converts it into hits per band at radius
%
% Input Arg:
% tab1 = zform table from pdist or pdist2
% dim = number of dimensions in euclidian space, if cross channel enter 2
% dr = edges of bands used in the radial density calculation
%
% Output Arg:
% a1 = hits within each band of the radial density function
%
function a1 = LAMA_CSR_distrange(tab1,dim,dr)
    if dim == 1
        sortedtab= sort(tab1,1);
    elseif dim == 2
        sizeans = size(tab1);
        [~,I] = max(sizeans);
        if I == 1
            sortedtab = sort(tab1,2)';
        elseif I == 2
            sortedtab = sort(tab1,1);
        end
    end
    a1 = zeros(3,width(dr)-1);
    for i = 1:(length(dr)-1)
       hits = tab1 >= dr(i)  ;
       hits = tab1(hits) < dr(i+1);
       a1(1,i) = sum(hits,'all');
       a1(2,i) = a1(1,i)/(size(tab1,1)*size(tab1,2));
       a1(3,i) = sum(tab1 < dr(i+1),'all')/(size(tab1,1)*size(tab1,2));
    end
end