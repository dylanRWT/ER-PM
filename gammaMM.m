function output = gammaMM(x,k,phase,high,prob,alpha,beta)
    %%
    %k = 2;
    %phase = [0,100];
    %high = [1,1];
    %prob=[0.6,0.4];
    %alpha = [10,10];
    %beta = [20,20];
    %output= NaN(length(x),k);
    for j = 1:k
        x1 = x-phase(1,j);
        z = pdf("Gamma",x1,alpha(1,j),beta(1,j));
        output(:,j) = (high(1,j)*z*prob(1,j))';
    end
    %%
end