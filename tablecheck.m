function [cleanArray,distNum] = tablecheck(maybetable,upperlimit)
    distNum = width(maybetable);
    if istable(maybetable)
        maybetable = table2array(maybetable);
    end
    % open each column as a data set and cull values over 2000 nm
    cleanArray = zeros(height(maybetable),distNum);
    for i = 1:distNum
        temp = maybetable(~isnan(maybetable(:,i)),i);
        temp = temp(temp <= upperlimit,1);
        d = height(maybetable) - height(temp);
        cleanArray(:,i) = [temp;zeros(d,1)];
    end
end