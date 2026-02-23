function corrected = unpadarray(padded,ogSize)
% unpadarray will take the padded dataset and crop it back to ogSize
% traditionally, unpadarray undos the kernel bleed at the edges of a
% dataset when you convolve.
%
% fchang@fas.harvard.edu

Bstart=ceil((size(padded)-ogSize)/2)+1;
Bend=Bstart+ogSize-1;
numDims = ndims(padded);
idx = cell(numDims,1);
for i = 1:numDims
   idx{i} = Bstart(i):Bend(i); 
end
corrected = padded(idx{:});
end

