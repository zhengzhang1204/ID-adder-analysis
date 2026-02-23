function [canvas,SList]=drawSptMapFLvImgFL_Nintv2(canvas,locList,Map,SList,PHMask,num)%in: ini(i.e., SList), cent(TS-scaled), map
% This is for Case 'adding spot' in ImgFL_N. Need to add a new line at the
% rear of existing SList, as well as a new nan line in map.

curLocList=locList(end-num+1:end,:);

for ix=1:num
%1. add a new line to the SList:
newCoor=curLocList(ix,:);
newSListRow=nan(1,9);

%Format of SList: [Frame, Idx, PercentageInCell, oriN, cellLength, SpotY(from image roof),SpotX(local tile), CoorInCell_Y]


Frame=floor(newCoor(1)/32)+1;

curMask=PHMask(:,(Frame-1)*32+1:Frame*32);
CC=bwconncomp(curMask,4);
L=regionprops(CC,'Centroid','PixelIdxList','PixelList');
[~,idxTmp]=sort(cellfun(@(x) x(2), {L.Centroid}),'descend');
X={L.PixelIdxList};
X=X(idxTmp);
Z={L.PixelList};
Z=Z(idxTmp);
newInd=sub2ind([256,32],newCoor(2),newCoor(1)-(Frame-1)*32);
Idx=find(cell2mat(cellfun(@(x) ismember(newInd, x), X,'UniformOutput',false)));
if isempty(Idx)
    disp('Added point is out of known cell regions.')
    return
end
newSListRow(1)=Frame;
newSListRow(2)=Idx;
curPixelList=Z{Idx};
newSListRow(5)=max(curPixelList(:,2))-min(curPixelList(:,2));
newSListRow(8)=max(curPixelList(:,2))-newCoor(2);
newSListRow(3)=newSListRow(8)/newSListRow(5);
SList=[SList;newSListRow];

%2. preparing out image with the added line.
canvas = insertMarker(canvas,curLocList(ix,:),'x-mark','color','green');
canvas = insertText(canvas,curLocList(ix,:), num2str(size(Map,1)),'BoxOpacity',0,'FontSize',9,'TextColor','green');
end
% out=cat(3,RJ+Contour,RJ+2*canvas(:,:,2),RJ);%with red contours and x-cross and spot id
% canvLineRawImg=canvas(:,:,1);%[uint8] black canvas with white lines only
end
