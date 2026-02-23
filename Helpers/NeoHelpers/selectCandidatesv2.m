function candidates = selectCandidatesv2(isPicOut,estimated,coe,varargin)
%SELECTCANDIDATES find candidate regions among the processed output of the
% data from findSpotsStage1
%
% detected:     output struct from findSpotsStage1
% candidates:   output data structure organized as ...
%
% [notes] - this function output candidates as connected components that
%           does not have to be rectangular bounding boxes.  this
%           accomodates the fitting of spots to complex shapes.
%         - this function filters the LLRATIO > thresh and A parameter > 0,
%           then smooths it to do hdome local maxima search.
%         - then for those connected objects smaller than some value gives
%
%
% [param cascade] -> simpleThresholdDetection


%--parameters--------------------------------------------------------------
% you can just use simple 'threshold' or 'hdome'
params.selectField          = 'LLRatio';       % selects which field to select on
params.strategy             = 'threshold';     % {'hdome','threshold','otsu'}
%==universal parameters====================================================
params.Athreshold           = 0;               % select regions where A > Athreshold
params.clearBorder          = true;            % clear border on xy perimeter
params.minVol               = 6;%was9              % min volume of feature. ZZ: important for low SNR images that contains a great number of small noise points.
params.imposeMinSize        = true;            % 
params.useMask              = [];              % use permissive mask
%==hdome specific parameters===============================================
params.hdomeH               = 1e3;
params.thresholdHDome       = 'threshold';          %{'otsu',thresholdValue}
%==threshold specific parameters===========================================
params.fieldThresh         = [];
%--------------------------------------------------------------------------
params = updateParams(params,varargin);

%% universal computation - preprocessing
if iscell(estimated.A1)
%     smoothField = estimated.convFunc(estimated.(params.selectField),estimated.spotKern{1});
    smoothField = estimated.(params.selectField);
    Athresholded = cellfunNonUniformOutput(@(x) x<params.Athreshold,estimated.A1);
    Athresholded = ~multiCellContents(Athresholded);
    
    if iscell(estimated.spotKern{1})
        sizeKern = cellfun(@(x) numel(x),estimated.spotKern{1});
    else
        sizeKern = size(estimated.spotKern{1});
    end
    
else %normally goes to here
    smoothField = estimated.convFunc(estimated.(params.selectField),estimated.spotKern);
    Athresholded = estimated.A1 > params.Athreshold;
    sizeKern = size(estimated.spotKern);
end


%% user specified computation
switch params.strategy
    case 'hdome'
        selectedRegions =  hdome(smoothField,params.hdomeH);
        if strcmp(params.thresholdHDome,'otsu')
            params.thresholdHDome = multithresh(selectedRegions(:));
        end
        selectedRegions = selectedRegions > params.thresholdHDome;
    case 'threshold'
        if isempty(params.fieldThresh)
            [params.fieldThresh, ~, ~] = threshold(multithresh(smoothField(:)), max(smoothField(:)), maxintensityproj(smoothField,3));
        end
        selectedRegions = smoothField > params.fieldThresh*coe;
    case 'otsu'
        thresh = multithresh(smoothField(:),1);
        selectedRegions = smoothField > thresh;
    otherwise
        error('unrecognized strategy');
end
%AsPicOut(selectedRegions, 76, 1,isPicOut);%ZZ added for debugging.---------------------1
%% universal computations - postprocessing
L = Athresholded.*selectedRegions;
%AsPicOut(L, 79, 1,isPicOut);%ZZ added for debugging.---------------------2
% use perssive mask if provided

if ~isempty(params.useMask)%<-------------------------change back to visible.
    L = maskData(L,params.useMask);   
end
%AsPicOut(L, 85, 1,isPicOut);%ZZ added for debugging.---------------------3
% remove small and not circular objects
L=ZZ_filterL(L,params.minVol);% was L = bwareaopen(L,params.minVol);
%AsPicOut(L, 88, 1,isPicOut);%ZZ added for debugging.---------------------4
if params.clearBorder
    L = clearXYBorder(L);
end
L92=L;
%AsPicOut(L, 92, 1,isPicOut);%ZZ added for debugging.---------------------5
if params.imposeMinSize
    L = L > 0;
    [~,~,seeds] = breakApartMasks(smoothField,L,varargin{:});
    seeds = seeds | Skeleton3D(L);
    minBBoxMask = imdilate(seeds,strel(ones(sizeKern(:)')));
    L = L | minBBoxMask;
end

[L,~,seeds] = breakApartMasks(smoothField,L>0,varargin{:});
%AsPicOut(L, 102, 1,isPicOut);%ZZ added for debugging.---------------------6
L = bwareaopen(L,params.minVol,4);%ZZ added 4 on 2021.7.20
%L = bwlabeln(L>0);
%AsPicOut(L, 105, 1,isPicOut);%ZZ added for debugging.---------------------7
stats = regionprops(L,estimated.(params.selectField),'PixelList','SubarrayIdx','PixelIdxList','BoundingBox','MaxIntensity','MinIntensity');
% remove any bbox that only has 1 dimension, this is not a spot
% remove any bbox that has the min = max values (this means the region is
% dead)
for ii = 1:numel(stats)
    currSize = genSizeFromBBox(stats(ii).BoundingBox);
    if any(currSize == 1)
         L(stats(ii).PixelIdxList) = 0;
    end
    if stats(ii).MinIntensity == stats(ii).MaxIntensity
        L(stats(ii).PixelIdxList) = 0;
    end
end
%L = bwareaopen(L>0,params.minVol,4);%ZZ added 4 on 2021.7.20
%AsPicOut(L, 120, 1,isPicOut);%ZZ added for debugging.---------------------8
L = L92 & L;
stats = regionprops(L,'PixelList','PixelIdxList','Centroid');
% need to have minimum volume
%candidates.L        = L;
% candidates.BWmask   = BWmask;
candidates.stats    = stats;
% candidates.seeds    = seeds;
% candidates.smoothField = smoothField;
candidates.params   = params;

