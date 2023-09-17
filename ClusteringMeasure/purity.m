function [Purity] = purity(predY,Y)
predLidx = unique(predY); 
pred_classnum = length(predLidx);
correnum = 0;
for ci = 1:pred_classnum
    incluster = double(Y(find(predY == predLidx(ci))));
    inclunub = hist(incluster, 1:max(incluster)); 
    if isempty(inclunub) inclunub=0;end;
    correnum = correnum + max(inclunub);
end;
Purity = correnum/length(predY);
end

