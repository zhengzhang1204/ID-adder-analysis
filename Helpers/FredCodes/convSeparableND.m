function dataTempConv = convSeparableND(data,separatedTemplate)
%CONVSEPARABLEFFTND executes separable convolution given that the template
%is in {templateX,templateY,templateZ,...} with each template oriented in
%the correct dimension.
%
% - unpads output to size(data)
% 
% fchang@fas.harvard.edu

sizeData=size(data);
numDims=length(sizeData);
dataTempConv = data;
for i = 1:numDims
    dataTempConv = convn(separatedTemplate{i},dataTempConv);
end

dataTempConv = unpadarray(dataTempConv,size(data));

