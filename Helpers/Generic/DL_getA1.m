function [out,estimated,me]=DL_getA1(data,spotKern,cameraVariance)
persistent spotKernSaved;
persistent cameraVarianceSaved;
persistent invVar;
persistent k1;
persistent k3;
persistent k5;

if iscell(spotKern)
    % convolution is separable
    convFunc = @convSeparableND;
    sqSpotKern = cellfunNonUniformOutput(@(x) x.^2,spotKern);
    onesSizeSpotKern = genSeparableOnes(cellfun(@(x) numel(x),spotKern));
else
    % otherwise just do fft
    convFunc = @convFFTND;
    sqSpotKern = spotKern.^2;
    onesSizeSpotKern = ones(size(spotKern));
end

% if spotKern || cameraVariance changes, update transformation matrix
if ~isequal(spotKernSaved,spotKern) || ~isequal(cameraVarianceSaved,cameraVariance)
    invVar = 1./cameraVariance;
    k1 = convFunc(invVar,spotKern);
    k3 = convFunc(invVar,sqSpotKern);
    k5 = convFunc(invVar,onesSizeSpotKern);
    spotKernSaved = spotKern;
    cameraVarianceSaved = cameraVariance;
end

me=[median(data,'all'),mean(data,'all')];

dataNormed = data.*invVar;
k2 = convFunc(dataNormed,spotKern);
k4 = convFunc(dataNormed,onesSizeSpotKern);
k6 = convFunc(dataNormed.*data,onesSizeSpotKern);

Normalization = k1.^2 - k5.*k3;

% parameters given A*spotKern + B, model of 1 spot
%A0          = k2./ k3;
A1          = (k1.*k4 - k5.*k2 ) ./ Normalization;
B1          = (k1.*k2 - k3.*k4)  ./ Normalization;
LL1         = -((B1.^2).*k5 + A1.*(2*B1.*k1 - 2*k2 + A1.*k3) - 2*B1.*k4 + k6);
% parmeters given B only, model of 0 spot
B0          = k4./k5;
LL0         = -((B0.^2).*k5 - 2*B0.*k4 + k6);

%estimated.A0         = unpadarray(A0,size(data));
estimated.A1         = unpadarray(A1,size(data));
%estimated.B1         = unpadarray(B1,size(data));
%estimated.B0         = unpadarray(B0,size(data));
estimated.LLRatio    = unpadarray(LL1-LL0,size(data));
estimated.spotKern   = spotKern;
out=estimated.A1;
% historical note: previous version of my code outputed A as Abefore noted
% below, which is incorrect.
% Abefore = (k2-k1.*k4./k5)./sqrt(k3./k5 - (k1./k5).^2);

estimated.convFunc   = convFunc;
end