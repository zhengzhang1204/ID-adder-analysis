function [daGauss,separateComponents] = ndGauss(sigmaSQVector,sizeVector,varargin)
%MAKE3DGAUSS generates a nD gaussian parameterized by:
% sigmaVector = [sigmasq_xx,sigmasq_y,sigmasq_z,....]
% sizeVector  = [numPixels_x,numPixels_y,numPixels_z,...]
% also a muVector can be defined as a last argument that defines the offset
% e.g. stuff = ndGauss([sigmas],[sizes],[offsets]);
% daGauss is the nD numeric array that contains the nD gaussian
% separateComponents is the separable components that can be used for
% separable convolution
% 
% [notes] - peak is set to 1.  if the dataset is even, then the invisible
% peak that is in between pixels is 1
%
%
% fchang@fas.harvard.edu

if nargin == 2
    muVector = zeros(size(sigmaSQVector));
else
    muVector = varargin{1};
end

dims = numel(sigmaSQVector);
separateComponents = cell(dims,1);
for i = 1:dims
    defineDomain = sizeVector(i)/2 - 0.5;
    currDomain = -defineDomain:defineDomain;
    currGauss = normpdf(currDomain,muVector(i),sqrt(sigmaSQVector(i)));
    % normalize each component, which normalizes total components
     currGauss = currGauss*sqrt(2*pi) .* sqrt(sigmaSQVector(i));
    reshapeVec = ones(dims,1);
    reshapeVec(i) = numel(currGauss);
    if isscalar(reshapeVec)
        
    else
        currGauss = reshape(currGauss,reshapeVec');
    end
    
    separateComponents{i} = currGauss;
    if i == 1
        daGauss = currGauss;
    else
        daGauss = convn(currGauss,daGauss);
    end
end


