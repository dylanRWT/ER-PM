function outputarray = tableconverter(x,y)
%%
    outputarray = [];
    for i = 1:length(x)
        if round(y(i)) > 0
        outputarray = [outputarray,repelem(x(i),round(y(i)))];
        end
    end
%%
end
    