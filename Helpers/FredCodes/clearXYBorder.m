function [cleared,touchingXY] = clearXYBorder(BWdata)
%CLEARXYBORDER clears positives values from the xy border of the data

xyBorder = true(size(BWdata));
coreSelect = cell(ndims(BWdata),1);
for ii = 1:ndims(BWdata)
    if ii > 2
        coreSelect{ii} = 1:size(BWdata,ii); 
    else
        coreSelect{ii} = 2:size(BWdata,ii)-1; 
    end
   
end
xyBorder(coreSelect{:}) = 0;

touchingXY = imreconstruct(xyBorder,BWdata);
cleared = BWdata;
cleared(touchingXY>0) = 0;
end

