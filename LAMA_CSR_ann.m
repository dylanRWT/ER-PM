% LAMA_CSR_ann
%
% Created by Dylan Greening for use in RWTurner Lab UCALGARY
%
% A function used in LAMA_CSR.m -> LAMA_CSR3.m
%
% Calculates average nearest neighbour from nearest neighbour arrays
%
%%
function ann = LAMA_CSR_ann(nnd)
    ann = zeros(size(nnd,1),size(nnd,2));
    for i = 1:height(nnd)
        for j = 1:width(nnd)
            ann(i,j) = mean(nnd{i,j});
        end
    end
end