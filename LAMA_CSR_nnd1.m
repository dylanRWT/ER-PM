%LAMA_CSR_nnd1
function nnd = LAMA_CSR_nnd1(pd_in,pooled)
nnd = {};
if pooled ~= 1
    for i = 1:height(pd_in)
        for j = 1:width(pd_in)
            nnd{i,j} = pd_in{i,j}(1,:)';
        end
    end
elseif pooled == 1
    nnd = [];
    for i = 1:height(pd_in)
        for j = 1:width(pd_in)
        temp = pd_in{i,j}(1,:)';
        nnd = [nnd;temp];
        end
    end
end