function [cIn]=DL_RefineCandv2(cIn, FLraw, Msk)
%1. create mask 
mskExSpots=true(size(Msk));
uit=zeros(length(cIn.stats),4);
for iCand=1:length(cIn.stats)
    mskExSpots(cIn.stats(iCand).PixelIdxList)=false;
end

%2. calculate background value
bkRaw_ExSp=median(nonzeros(double(FLraw).*mskExSpots));% median of non-cell background

%3. record mean/bk, median/bk, Q95/bk, area, center.
for iCand2=1:length(cIn.stats)
    uit(iCand2,1)=mean(FLraw(cIn.stats(iCand2).PixelIdxList),'all')./bkRaw_ExSp;
    uit(iCand2,2)=median(FLraw(cIn.stats(iCand2).PixelIdxList),'all')./bkRaw_ExSp;
    uit(iCand2,3)=quantile(FLraw(cIn.stats(iCand2).PixelIdxList),0.98)./bkRaw_ExSp;
    uit(iCand2,4)=length(cIn.stats(iCand2).PixelIdxList);
    uit(iCand2,5:6)=round(median(cIn.stats(iCand2).PixelList,1));
end
cIn.IntInfo=uit;
end


