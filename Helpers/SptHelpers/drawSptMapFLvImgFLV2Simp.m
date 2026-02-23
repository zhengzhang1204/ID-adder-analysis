
function [out]=drawSptMapFLvImgFLV2Simp(locList,RJ,PHMask)
% Preparation and initation:
Contour=uint8(imdilate(PHMask,strel('diamond',2))-PHMask)*255*0.4;
canvas=uint8(zeros(size(PHMask)));

% remove bad entries for spots.
del=find(isnan(locList(:,1)));
if ~isempty(del)
    locList(del,:)=[];
end

% draw marker:
canvas = insertMarker(canvas,locList,'x-mark','color','green');

out=cat(3,RJ+Contour,RJ+2*canvas(:,:,2),RJ);
end
