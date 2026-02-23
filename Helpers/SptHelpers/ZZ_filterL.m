function [outL]=ZZ_filterL(L,minVol)
L=bwareaopen(L,minVol);
bwL=bwlabel(L,4);
stats=regionprops(bwL,'Area','Eccentricity','PixelIdxList');
for i=1:length(stats)
    if stats(i).Area<15 && stats(i).Eccentricity>0.7
        bwL(stats(i).PixelIdxList)=0;
    end
end
outL=logical(bwL);
end