function [ masked ] = maskData( data,mask )
%MASKDATA Summary of this function goes here
%   Detailed explanation goes here

sizeData = size(data);
sizeMask = zeros(ndims(data),1);
sizeMask(1:ndims(mask)) = size(mask);

extrusionFactor =(sizeData(:)-sizeMask(:));
extrusionFactor = extrusionFactor-1;
extrusionFactor(extrusionFactor<0) = 0;
extrudedMask = padarray(mask,extrusionFactor,'pre','replicate');
masked = data.*extrudedMask;


end

