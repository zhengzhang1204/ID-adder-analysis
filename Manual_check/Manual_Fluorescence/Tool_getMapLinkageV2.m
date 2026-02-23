function [uit]=Tool_getMapLinkageV2(idxM,idxD,PHmap)
% uit: 
% [0]: idxM == idxD;
% [1xN double]: [idxM, linker(s)];
% []: no linkage between idxM and idxD.
uit=[];
if idxM==idxD
    uit=0;
    return;
end
while idxD~=idxM
    if isnan(idxD)
        uit=[];%No link. Output is empty.
        return;
    end
    target=PHmap(idxD,1);
    if target~=0
        uit=cat(2,target,uit);
        idxD=target;
    else
        uit=[];
        return;
    end
end
% if ~isempty(uit)
%     uit(1)=[];
% end
end