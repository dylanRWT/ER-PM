% LAMA_CSR_err
%
% Created by Dylan Greening for use in RWTurner lab UCalgary 
%
% A small function which takes data as a array (so that entries (columns)
% can be different lengths) calculates the mean of each column along with
% the SEM and outputs these values for quick plotting and reference
%
function [mu,errhigh,errlow,sem] = LAMA_CSR_err(data)

    errhigh = [];
    errlow = [];
    mu = [];
    sem = [];
    for col = 1:width(data)
        mu(1,col) = mean(data{1,col});
        n = height(data{1,col});
        sem(1,col) = std(data{1,col})./sqrt(n);
        errhigh(1,col) = mu(1,col)+sem(1,col);
        errlow(1,col) = mu(1,col)-sem(1,col);   
    end

end