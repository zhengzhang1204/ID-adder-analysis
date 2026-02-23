function Ld = doWaterShedSplitting(dataND,seeds,mask)
%DOWATERSHEDSPLITTING 

D2 = imimposemin(-dataND,seeds);
Ld = watershed(D2);
Ld(mask==0) = 0;