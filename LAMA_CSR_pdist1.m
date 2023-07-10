%% LAMA_CSR_pdist1
%
% Created by Dylan Greening for use in RWTurner lab UCalgary 
%
% A function for use in LAMA_CSR.m -> LAMA_CSR3.
% 
% Used to define euclidian distance in from pdist()
%
function pd_out = LAMA_CSR_pdist1(cent)
    pd_out = {};
    for i = 1:height(cent)
        for j = 1:width(cent)
            temp = sort(squareform(pdist(cent{i,j})),1);
            pd_out{i,j} = temp(2:end,:);
        end
    end
end