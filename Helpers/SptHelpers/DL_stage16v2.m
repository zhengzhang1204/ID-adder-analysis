function [Cand,Est,Msk]=DL_stage16v2(Thr,Est,Mask,TSstack,minArea,COE,isPicOut,Cand)%,cropInf,imgSize,YcropF)
%addons={'cellMaskVariable','genMaskWOtsu1','cellMasks',[],'selectField','LLRatio','doProcParallel',true};
nSchn=size(TSstack,3);
nTP=length(Mask);

%1. Generate tilized stack (TS) for masks:
Msk=false(256,32*nTP,nSchn);
for iTP=2:length(Mask)
    for iS=1:nSchn
        Msk(:,(iTP-1)*32+1:iTP*32,iS)=Mask{1,iTP}{iS,1};
    end
end

%2. Process each side channel:
if isempty(Cand); Cand=cell(nSchn,1); end
for iSchn=1:nSchn%PARFOR each side channels<--------------change this back to 1
% Register empty side channels.
    if isempty(find(Msk(:,:,iSchn),1))
        Est{1,iSchn}=[];
        continue;
    end
% Retrieve candidates.
    if isempty(Est{1,iSchn})
        continue;
    end
    CandTemp = selectCandidatesv2(isPicOut,Est{1,iSchn} ,COE,'useMask',Msk(:,:,iSchn),...
    'fieldThresh',Thr(iSchn,1),'cellMaskVariable','genMaskWOtsu1','cellMasks',[],...
    'selectField','LLRatio','doProcParallel',true);

% Modules to refine candidates:
% Module-1 remove small objects and reduce size of Est:
    CandTemp = remSmallCand(CandTemp,minArea);
    if length(CandTemp.stats)<10
        Est{1,iSchn}=[];
        continue;
    end
% Module-2 Retrieve FL information from FLraw images in struct Cand 
    [CandTemp]= DL_RefineCandv2(CandTemp,TSstack(:,:,iSchn),Msk(:,:,iSchn));
    Cand{iSchn,1}=CandTemp;
end
end
function [uit] = remSmallCand(in,thres)
delSel=[];
for i=1:length(in.stats)
    if length(in.stats(i).PixelIdxList)<thres
        delSel=cat(1,delSel,i);
    end
end
in.stats(delSel)=[];
%in = rmfield(in,{'L','seeds','smoothField'});
uit=in;
end