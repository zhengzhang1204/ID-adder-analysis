function [uit]=Tool_getMapLinkage(idxM,idxD,PHmap)
uit=[];
while idxD~=idxM
    if isnan(idxD)
        uit=0;%meaning the mother and daughter has no link.
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
if ~isempty(uit)
    uit(1)=[];
end
end