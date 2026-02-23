
function [out]=ImgFLN_removeOldSpot(RJ,SpotCent,PHMask)
Contour=uint8(imdilate(PHMask,strel('diamond',2))-PHMask)*255*0.4;
canvas=uint8(zeros(size(PHMask)));
nFrm=size(RJ,2)/32;

for iLin=1:size(SpotCent,1)
% draw markers and spot row ID:
    if ~isnan(SpotCent(iLin,1))
        canvas = insertMarker(canvas,SpotCent(iLin,:),'x-mark','color','green');
        canvas = insertText(canvas,SpotCent(iLin,:), num2str(iLin),'BoxOpacity',0,'FontSize',9,'TextColor','green');
    end
end

%draw frame numbers:
TxtPos=[((0:1:nFrm-1)*32+6)',235*ones(nFrm,1)];
txt=num2cell((1:1:nFrm)');
txt=cellfun(@(x) num2str(x),txt,'UniformOutput',false);
canvas = insertText(canvas,TxtPos,txt,'BoxOpacity',0,'FontSize',16,'TextColor','white');
out=cat(3,RJ+Contour,RJ+2*canvas(:,:,2),RJ);
end
